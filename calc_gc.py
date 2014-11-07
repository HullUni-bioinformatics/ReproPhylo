def calc_gc(sequence):
    sequence = sequence.upper()                    # make all chars uppercase
    n = sequence.count('T') + sequence.count('A')  # count only A, T,
    m = sequence.count('G') + sequence.count('C')  # C, and G -- nothing else (no Ns, Rs, Ws, etc.)
    if n + m == 0:
        return 0.                                  # avoid divide-by-zero
    return float(m) / float(n + m)

def test_1():
    result = round(calc_gc('ATGGCAT'), 2)
    print 'hello, this is a test; the value of result is', result
    assert result == 0.43
    
def test_2(): # test handling N
    result = round(calc_gc('NATGC'), 2)
    assert result == 0.5, result
    
def test_3(): # test handling lowercase
    result = round(calc_gc('natgc'), 2)
    assert result == 0.5, result