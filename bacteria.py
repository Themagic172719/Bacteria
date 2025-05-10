from fastaReader import fastaReader
import random
import numpy as np
import copy
from evaluadorBlosum import evaluadorBlosum

class bacteria:
    
    def __init__(self, path):
        self.matrix = fastaReader(path)
        self.blosumScore = 0
        self.fitness = 0
        self.interaction = 0
        self.NFE = 0

    def showGenome(self):
        for seq in self.matrix.seqs:
            print(seq)

    def clonar(self, path):
        newBacteria = bacteria(path)
        newBacteria.matrix.seqs = np.array(self.matrix.seqs, copy=True)
        return newBacteria

    def tumboNado(self, numGaps):
        self.cuadra()
        matrixCopy = np.array(self.matrix.seqs, copy=True).tolist()
        for _ in range(random.randint(0, numGaps)):
            seqnum = random.randint(0, len(matrixCopy) - 1)
            pos = random.randint(0, len(matrixCopy[0]))
            matrixCopy[seqnum] = matrixCopy[seqnum][:pos] + "-" + matrixCopy[seqnum][pos:]
        self.matrix.seqs = np.array(matrixCopy)
        self.cuadra()
        self.limpiaColumnas()

    def cuadra(self):
        maxLen = max(len(s) for s in self.matrix.seqs)
        self.matrix.seqs = np.array([s.ljust(maxLen, '-') for s in self.matrix.seqs])

    def gapColumn(self, col):
        return all(seq[col] == "-" for seq in self.matrix.seqs)

    def limpiaColumnas(self):
        i = 0
        while i < len(self.matrix.seqs[0]):
            if self.gapColumn(i):
                self.deleteCulmn(i)
            else:
                i += 1

    def deleteCulmn(self, pos):
        for i in range(len(self.matrix.seqs)):
            self.matrix.seqs[i] = self.matrix.seqs[i][:pos] + self.matrix.seqs[i][pos + 1:]

    def getColumn(self, col):
        return [seq[col] for seq in self.matrix.seqs]

    def autoEvalua(self):   
        evaluador = evaluadorBlosum()
        score = 0
        for i in range(len(self.matrix.seqs[0])):
            column = self.getColumn(i)
            gapCount = column.count("-")
            column = [x for x in column if x != "-"]
            pares = self.obtener_pares_unicos(column)
            for par in pares:
                score += evaluador.getScore(par[0], par[1])
            score -= gapCount * 2
        self.blosumScore = score
        self.NFE += 1

    def obtener_pares_unicos(self, columna):
        return list({tuple(sorted([columna[i], columna[j]])) for i in range(len(columna)) for j in range(i + 1, len(columna))})
