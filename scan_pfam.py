import sys
import requests
import xml.etree.ElementTree as ET
import io
import csv
from Bio import SeqIO
import xml

def send_request(request_string):
    """Returns uuids (id of the query)"""
    url = 'https://www.ebi.ac.uk/Tools/hmmer/search/hmmscan'
    body = {'seq': request_string, 'hmmdb': 'pfam'}
    headers = {'Accept': 'text/xml', 'Expect': None}
    r = requests.post(url, body, headers=headers)
    try:
        return ET.fromstring(r.text).find('data').get('uuid')
    except AttributeError:
        print(f'UUID not found: {request_string}')
        return None
    except xml.etree.ElementTree.ParseError:
        print(f'UUID not found: {request_string}')
        return None

def receive_tsv(uuid):
    url = f'https://www.ebi.ac.uk/Tools/hmmer/download/{uuid}/score'
    body = {'format': 'tsv'}
    r = requests.post(url, body)
    return r.text

def extract_domains_family_ids(tsv):
    reader = csv.DictReader(io.StringIO(tsv), dialect='excel-tab')
    ids = []
    try:
        for row in reader:
            ids.append(row['Family id'])
        return ids
    except AttributeError:
        print(f'Family id not found', file=sys.stderr)
        return []

def make_array(occurences_dict):
    def flatten(l):
        return [item for sublist in l for item in sublist]

    all_domains = list(set(flatten(occurences_dict.values())))
    yield ['ID'] + all_domains
    for key, domains in occurences_dict.items():
        yield ([key] + [1 if gdomain in domains else 0 for gdomain in all_domains])

def write_csv(filename, generator):
    with open(filename, 'w') as f:
        filewriter = csv.writer(f, delimiter=',')
        for line in generator:
            filewriter.writerow(line)

"""MAIN"""

if len(sys.argv) == 3:
    filename_in = sys.argv[1]
    filename_out = sys.argv[2]
else:
    print("Usage: python3 scan_pfam.py input_filename output_filename", file=sys.stderr)
    exit(0)

print(f'Opening file to read: {filename_in}')
sequences = list(SeqIO.parse(filename_in, "fasta"))
domains = dict()

print(f'Starting searching for the domains')
slen = len(sequences)
i = 1
errors = 0
for s in sequences:
    print(f'Sending hmmscan request for: {s.id}: {i}/{slen}')
    i = i + 1
    uuid = send_request(str(s.seq))
    print(f'Received UUID: {uuid}')
    if uuid == None:
        print(f'Skipped {s.id}')
        errors = errors + 1
        continue
    print(f'Downloading tsv, uuid: {uuid}')
    tsv = receive_tsv(uuid)
    print(f'Parsing tsv')
    d = extract_domains_family_ids(tsv)
    domains[s.id] = d

print(f'Saving results to csv file: {filename_out}')
print(f'{errors} queries failed, probably due to malformed xml')
write_csv(filename_out, make_array(domains))
