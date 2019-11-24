import requests
import re
import time
import sys
import xml.etree.ElementTree as ET
from Bio import SeqIO
import io


def send_request(request_string):
    """Sends request to blastp, returns RID."""
    url = 'https://blast.ncbi.nlm.nih.gov/blast/Blast.cgi'
    body = {'CMD':'Put',
            'QUERY': request_string,
            'DATABASE': 'nr',
            'PROGRAM': 'blastp',
            'THRESHOLD': '21',
            'WORD_SIZE': '6',
            'COMPOSITION_BASED_STATISTICS': '2',
            'ALIGNMENTS': '100',
            'MATRIX': 'BLOSUM62',
            'GAPCOSTS': '11 1',
            'EXPECT': '10',
            'HITLIST': '100',
            }
    r = requests.post(url, body)
    pattern = re.compile(r'QBlastInfoBegin.*RID = (.*?)\n.*QBlastInfoEnd', re.DOTALL)

    return pattern.findall(r.text)[0]


def get_status(rid):
    """Gives status of RID. It's either: WAITING, READY or UNKNOWN"""
    url = 'https://blast.ncbi.nlm.nih.gov/blast/Blast.cgi'
    r = requests.get(url, {'CMD': 'Get', 'RID': rid})
    pattern = re.compile(r'QBlastInfoBegin.*Status=(.*?)\n.*QBlastInfoEnd', re.DOTALL)

    return pattern.findall(r.text)[0]


def receive_xml(rid):
    status = get_status(rid)
    waiting = 0
    while status == 'WAITING':
        time.sleep(15)
        waiting += 15
        status = get_status(rid)
        print(f'Waiting for {waiting} seconds\r', end='', file=sys.stderr)
    assert status == 'READY'

    url = 'https://blast.ncbi.nlm.nih.gov/blast/Blast.cgi'
    r = requests.get(url, {'CMD': 'Get', 'RID': rid, 'FORMAT_TYPE': 'XML'})
    return r.text


def parse_xml(text):
    root = ET.fromstring(text)
    iterations = root.find('BlastOutput_iterations').findall('Iteration')
    ret = []

    for iteration in iterations:
        hits = iteration.find('Iteration_hits')

        for hit in hits.findall('Hit'):
            align_ok = float(hit.find('Hit_hsps').find('Hsp').find('Hsp_positive').text)
            align_len = float(hit.find('Hit_hsps').find('Hsp').find('Hsp_align-len').text)

            accession_number = hit.find('Hit_accession').text
            similarity = align_ok/align_len
            evalue = float(hit.find('Hit_hsps').find('Hsp').find('Hsp_evalue').text)
            ret.append((accession_number, similarity, evalue))
    return ret


def filter_proteins(proteins, e_similarity= 0.9, e_evalue=10**(-10)):
    ret = []
    for protein in proteins:
        _, similarity, evalue = protein
        if e_similarity <= similarity and evalue <= e_evalue:
            ret.append(protein)
    return ret


def get_protein_names(proteins):
    ret = []
    for accession_number, _, _ in proteins:
        ret.append(accession_number)
    return ret


def accession_numbers_to_SeqRecords(accession_numbers):
    """
    accession_numers: list of strings of accesions numbers
    of proteins (e.g. ['WP_054314356', 'EUJ37427'])
    returns list of SeqRecords

    Uses https://eutils.ncbi.nlm.nih.gov/entrez/eutils
    to get protein sequences from accesion numbers
    or None if it's not in the database.
    """

    uids = []
    url = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi'
    # using search engine to find proper uids (ids in the database)
    for an in accession_numbers:
        body = {'db': 'protein', 'term': an}
        try:
            uids.append(ET.fromstring(requests.get(url, body).text).find('IdList').find('Id').text)
        except AttributeError:
            print(f'Protein: {an} not found in the eutils database', file=sys.stderr)

    # retrieveing the sequences using the uids
    url = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi'
    body = {'db': 'protein', 'id': ','.join(uids), 'rettype': 'fasta'}

    r = requests.post(url, body)

    # using request result string as a file
    ret = list(SeqIO.parse(io.StringIO(r.text), "fasta"))
    if len(ret) != len(accession_numbers):
        print(f"accession_numbers_to_SeqRecords: not all proteins found, {len(ret)} != {len(accession_numbers)}", file=sys.stderr)
    return ret



"""MAIN"""

args_message = """Usage: python3 extend.py filename [similarity evalue],
       similarity in percentages,
       evalue in python scientific notation (eg. 8e-2 = 8 * 10^(-2))
       filename.fasta should exist in the directory

       for example:
       python3 extend.py input.fasta 90 10e-10
"""

if len(sys.argv) == 1:
    filename = "input-z2.fasta"
    similarity = 0.9
    evalue = 10.0**(-10)
elif len(sys.argv) == 4:
    filename = sys.argv[1]
    similarity = float(sys.argv[2])/100
    evalue = float(sys.argv[3])
else:
    print(args_message, file=sys.stderr)
    exit(0)

print(f'Opening file: {filename}')
with open(filename) as f:
    input_data = f.read()

print(f'Sending request to blastp: bytes: {len(input_data)}')
rid = send_request(input_data)

print(f'Waiting for response, RID: {rid}')
blastp_xml = receive_xml(rid)

print(f'Received response, RID: {rid}')
print(f'Saving to backup file: blastp_backup.xml')
with open("blastp_backup.xml", "w") as f:
    f.write(blastp_xml)

print(f'Parsing xml file')
proteins = parse_xml(blastp_xml)

print(f'Filtering proteins: evalue: {evalue}, similarity: {similarity}')
proteins = filter_proteins(proteins, similarity, evalue)

print(f'Extracting accession numbers from filtered proteins')
protein_names = list(set(get_protein_names(proteins)))

print('Accessing database to find exact sequences from accession_numbers')
seq_records = accession_numbers_to_SeqRecords(protein_names)

print('Saving final results to blastp_result.fasta')
with open("blastp_result.fasta", "w") as f:
        SeqIO.write(seq_records, f, "fasta")
