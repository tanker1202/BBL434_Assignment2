import csv

def read_fasta(path):
    seq = ""
    with open(path) as f:
        for line in f:
            if not line.startswith(">"):
                seq += line.strip().upper()
    return seq


def read_scoring(csv_file):

    with open(csv_file) as f:
        reader = list(csv.reader(f))

    # first two rows = gap values
    gap_open = int(reader[0][1])
    gap_extend = int(reader[1][1])

    # third row = header
    bases = reader[2][1:]

    matrix = {}

    # matrix starts from row 3 onward
    for row in reader[3:]:
        rbase = row[0]
        for cbase, val in zip(bases, row[1:]):
            matrix[(rbase,cbase)] = int(val)

    return gap_open, gap_extend, matrix


def write_output(score, a1, a2):
    with open("Output.fa","w") as f:
        f.write(">Alignment_Sequence_1\n"+a1+"\n")
        f.write(">Alignment_Sequence_2\n"+a2+"\n")


def smith_waterman_affine(seq1, seq2, gap_open, gap_extend, matrix):

    n = len(seq1)
    m = len(seq2)

    M = [[0]*(m+1) for _ in range(n+1)]
    X = [[0]*(m+1) for _ in range(n+1)]
    Y = [[0]*(m+1) for _ in range(n+1)]

    ptrM = [[None]*(m+1) for _ in range(n+1)]
    ptrX = [[None]*(m+1) for _ in range(n+1)]
    ptrY = [[None]*(m+1) for _ in range(n+1)]

    max_score = 0
    max_pos = (0,0)
    max_matrix = "M"

    for i in range(1,n+1):
        for j in range(1,m+1):

            s = matrix[(seq1[i-1], seq2[j-1])]

            # gap in seq2
            open_gap = M[i-1][j] + gap_open
            extend_gap = X[i-1][j] + gap_extend
            X[i][j] = max(0, open_gap, extend_gap)

            if X[i][j] == 0:
                ptrX[i][j] = None
            elif X[i][j] == open_gap:
                ptrX[i][j] = ("M", i-1, j)
            else:
                ptrX[i][j] = ("X", i-1, j)

            # gap in seq1
            open_gap = M[i][j-1] + gap_open
            extend_gap = Y[i][j-1] + gap_extend
            Y[i][j] = max(0, open_gap, extend_gap)

            if Y[i][j] == 0:
                ptrY[i][j] = None
            elif Y[i][j] == open_gap:
                ptrY[i][j] = ("M", i, j-1)
            else:
                ptrY[i][j] = ("Y", i, j-1)

            # match matrix
            candidates = [
                (0,None),
                (M[i-1][j-1] + s, ("M", i-1, j-1)),
                (X[i-1][j-1] + s, ("X", i-1, j-1)),
                (Y[i-1][j-1] + s, ("Y", i-1, j-1))
            ]

            M[i][j], ptrM[i][j] = max(candidates, key=lambda x:x[0])

            # track max
            for name,val in [("M",M[i][j]),("X",X[i][j]),("Y",Y[i][j])]:
                if val > max_score:
                    max_score = val
                    max_pos = (i,j)
                    max_matrix = name


    aligned1 = []
    aligned2 = []

    i,j = max_pos
    matrix_name = max_matrix

    while True:

        if matrix_name == "M":
            if M[i][j]==0 or ptrM[i][j] is None:
                break
            prev = ptrM[i][j]
            aligned1.append(seq1[i-1])
            aligned2.append(seq2[j-1])

        elif matrix_name == "X":
            if X[i][j]==0 or ptrX[i][j] is None:
                break
            prev = ptrX[i][j]
            aligned1.append(seq1[i-1])
            aligned2.append("-")

        else:
            if Y[i][j]==0 or ptrY[i][j] is None:
                break
            prev = ptrY[i][j]
            aligned1.append("-")
            aligned2.append(seq2[j-1])

        matrix_name,i,j = prev

    return max_score,"".join(reversed(aligned1)),"".join(reversed(aligned2))


if __name__ == "__main__":

    seq1 = read_fasta("Input1.fa")
    seq2 = read_fasta("Input2.fa")

    gap_open, gap_extend, matrix = read_scoring("scoring.csv")

    score,a1,a2 = smith_waterman_affine(
        seq1,seq2,
        gap_open,gap_extend,
        matrix
    )

    print("\nBest Local Alignment\n")
    print("Score:",score)
    print(a1)
    print(a2)

    write_output(score,a1,a2)
