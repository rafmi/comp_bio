import csv
import pandas as pd
import math
import scipy.stats
import sys


def stirling_ln(n):
    '''calculates ~ln(n!)'''
    if n == 0:
        return float(0)
    if n < 0:
        print('Warning, n<0 in stirling')
        return float(0)
    if n < 1000:
        return math.log(math.factorial(n))
    n = float(n)
    a = math.log(2 * math.pi)/2
    b = math.log(n) * (n + 0.5)
    c = -n
    d = 1/(12 * n)
    e = -1/(360 * n**3)
    f = 1/(1260*n**5)
    e = -1/(1680*n**7)
    return a + b + c + d + e + f + e

def binomial_ln(a, b):
    if a < 0 or b < 0:
        print('Warning, a < 0 or b < 0 in binomial_ln')
    if b > a:
        print('warning b > a in binomial_ln')
        return 0
    if b == a or b == 0:
        return 0
    if a < 30:
        return math.log(math.factorial(a)/math.factorial(b)/math.factorial(a-b))
    x = min(b, a - b)
    y = max(b, a - b)
    ret = float(0)
    i = a
    while i > y:
        ret = ret + math.log(i)
        i = i - 1
    ret = ret - stirling_ln(x)

    return ret

def fisher_element1(i, n, K, N):
    a = binomial_ln(K, i)
    b = binomial_ln(N-K, n-i)
    c = binomial_ln(N, n)

    return math.exp(a-c+b)

def divide_factorials(up, down):
    '''First remembers every factor of every factorial in a dictionary
    and reduces equal factors
    then calculates factorials from up / factorials from down in safe order'''
    enc = {}
    for u in up:
        for i in range(1, u+1):
            enc[i] = enc.get(i, 0) + 1
    for d in down:
        for i in range(1, d+1):
            enc[i] = enc.get(i, 0) - 1

    times = []
    div = []

    for k, v in enc.items():
        if v > 0:
            times += [k] *  v
        elif v < 0:
            div += [k] * (-v)

    ret = 1.0
    times.sort()
    div.sort()
    while times != [] or div != []:

        if ret > 1:
            if div != []:
                ret = ret / div.pop()
            else:
                ret = ret * times.pop()
        else:
            if times != []:
                ret = ret * times.pop()
            else:
                ret = ret / div.pop()
    return ret

def fisher_element2(i, n, K, N):
    '''divides factorials into two groups: up(numerator), down(denumerator)'''
    up = (K, N-K, n, N-n)
    down = (i, K-i, n-i, N, N-K-n+i)
    return divide_factorials(up, down)

def fisher(fisher_element, i, n, K, N):
    return min(1.0, sum(fisher_element(j, n, K, N) for j in range(i, min(n, K, N)+1)))

def matches_from_file(filename):
    print(f'Analyzing {filename}')
    ret = dict()
    df = pd.read_csv(filename)
    length = len(df[df.columns[0]])
    for colname in df.columns[1:]:
        matches = sum(df[colname])
        print(f'{colname} occured in {matches} proteins in {filename}')
        ret[colname] = matches
    return ret, length

def merge_matches(matches1, matches2):
    ret = {}
    keys = set(matches1.keys()) | set(matches2.keys())
    for key in keys:
        ret[key] = (matches1.get(key, 0), matches2.get(key, 0))
    return ret

if len(sys.argv) != 3:
    message = "Usage: python3 fisher.py big_file.csv small_file.csv"
    print(message, file=sys.stderr)
    exit(0)

filename_big = sys.argv[1]
filename_small = sys.argv[2]

matches1, pcount1 = matches_from_file(filename_big)
print()
matches2, pcount2 = matches_from_file(filename_small)
matches = merge_matches(matches1, matches2)

print('\n\nFISHER TEST\n\n')
significant = []
for fname, (m1, m2) in matches.items():
    N = pcount1 + pcount2
    K = pcount1
    n = m1 + m2
    i = m1
    #table = [[m1, m2], [pcount1-m1, pcount2-m2]]
    table = [[m1, pcount1-m1], [m2, pcount2-m2]]

    f1 = fisher(fisher_element1, i, n, K, N)
    f2 = fisher(fisher_element2, i, n, K, N)
    (_, fsg) = scipy.stats.fisher_exact(table, 'greater')
    (_, fsd) = scipy.stats.fisher_exact(table)

    print(f'++++++++++++++++++++++++++++++++')
    print(f'CONTINGENCY TABLE for {fname}:\n')
    print(pd.DataFrame([
        ['', 'match', 'no match'],
        [filename_big, str(table[0][0]), str(table[0][1])],
        [filename_small, str(table[1][0]), str(table[1][1])]]).to_string(index=False, header=False))
    print()


    print(f'Fisher test implementation 1 (greater): {f1}')
    print(f'Fisher test implementation 2 (greater): {f2}')
    print('iterations:', min(n, K, N)+1 - i)
    print(f'Fisher test scipy greater: {fsg}')
    print(f'Fisher test scipy two-sided(default): {fsd}')
    if f2 <= 0.05:
        significant.append(fname)
        print(f'{fname} seems significant!!!\n')
    if fsd/fsg < 0.1:
        if (1-f2) <= 0.05:
            significant.append(fname)
        print(f'{fname} wrong assumption about sizes, is it {1-f2}?\n')
    else:
        print()

print('Some interesting domains:')
print('-------------------')
for sig in significant:
    f1m, f2m = matches[sig]
    print(f"Matches in {filename_big} for {sig}: {f1m}/{pcount1}")
    print(f"Matches in {filename_small} for {sig}: {f2m}/{pcount2}\n")
print('They seem to contradict the assumption about\neven distribution of protein domains in both files')
