with open('dataset.txt', 'r') as f:
    Text = f.readline()
    vals = f.readline().split()

Text = Text[:-1]

k = int(vals[0])
d = int(vals[1])

depth = k
FrequencyArray = []

kmers = ['' for i in range(4**k)]

def getLetter(num):
    return {
        0: 'A',
        1: 'C',
        2: 'G',
        3: 'T'
        }.get(num)

def getKmers(k, depth):
    if depth == 0:
        return
    iters = 4**(k-depth)
    count = 0
    for j in range(0, iters):
        for i in range(0, (4**depth)/4):
            kmers[count] += 'A'
            count += 1
        for i in range((4**depth)/4, 2*(4**depth)/4):
            kmers[count] += 'C'
            count += 1
        for i in range(2*(4**depth)/4, 3*(4**depth)/4):
            kmers[count] += 'G'
            count += 1
        for i in range(3*(4**depth)/4, 4**depth):
            kmers[count] += 'T'
            count += 1
    getKmers(k, depth-1)

def SymbolToNumber(symbol):
    if symbol == 'A':
        return 0
    if symbol == 'C':
        return 1
    if symbol == 'G':
        return 2
    if symbol == 'T':
        return 3
    
def PatternToNumber(pattern):
    if pattern == '':
        return 0
    symbol = pattern[-1]
    pattern = pattern[:-1]
    
    return 4*PatternToNumber(pattern) + SymbolToNumber(symbol)

def NumberToPattern(index):
    return kmers[index]

def HammingDistance(string1, string2):
    count = 0
    for i in range(0, len(string1)):
        if string1[i] != string2[i]:
            count += 1
    return count

def Suffix(Pattern):
    return Pattern[1:]

def Neighbors(Pattern, d):
    results = []
    if d == 0:
        return Pattern
    elif len(Pattern) == 1:
        results.append('A')
        results.append('C')
        results.append('G')
        results.append('T')
        return results
    Neighborhood = []
    SuffixNeighbors = Neighbors(Suffix(Pattern), d)
    for Text in SuffixNeighbors:
        if HammingDistance(Suffix(Pattern), Text) < d:
            Neighborhood.append('A' + Text)
            Neighborhood.append('C' + Text)
            Neighborhood.append('G' + Text)
            Neighborhood.append('T' + Text)
        else:
            Neighborhood.append(Pattern[0] + Text)
    return Neighborhood

def ComputeFrequencies(Text, k):
    for i in range(0, 4**k):
        FrequencyArray.append(0)
        
    for i in range(0, len(Text) - k):
        pattern = Text[i: i+k]
        nearbies = Neighbors(pattern, d)
        for kmer in nearbies:
            if HammingDistance(kmer, pattern) <= d:
                FrequencyArray[PatternToNumber(kmer)] += 1
        nearbies = Neighbors(ReverseComplement(pattern), d)
        for kmer in nearbies:
            if HammingDistance(kmer, ReverseComplement(pattern)) <= d:
                FrequencyArray[PatternToNumber(kmer)] += 1

    count = 0
    for item in FrequencyArray:
        if NumberToPattern(count) == ReverseComplement(NumberToPattern(count)):
            item = item/2
        count += 1


def MostFrequent():
    most = max(FrequencyArray)
    count = 0
    indeces = []
    for item in FrequencyArray:
        if item == most:
            indeces.append(count)
        count += 1
    return indeces

def ReverseComplement(String):
    Complement = ""
    for x in String:
        if x == 'A':
            Complement += 'T'
        elif x == 'C':
            Complement += 'G'
        elif x == 'G':
            Complement += 'C'
        elif x == 'T':
            Complement += 'A'

    return Complement[::-1]

if k <= 12:
    filename = str(k) + '.txt'
    with open(filename, 'r') as f:
        kmers = f.read().splitlines()
else:
    getKmers(k, k)
        
ComputeFrequencies(Text, k)
indeces = MostFrequent()
results = ""

for index in indeces:
    results += NumberToPattern(index) + " "

print(results)
