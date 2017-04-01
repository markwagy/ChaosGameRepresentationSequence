# Chaos Game Representation #

import matplotlib.pyplot as plt
import numpy as np
from PIL import Image
from Bio import SeqIO
import pickle

class CGR:

    def __init__(self, size=512, k=9):
        self.imgSize = size
        # k is the number of genes to consider at one time
        self.k = k

    def charPoint(self, char):
        if char == "A":
            return 0, 0
        if char == "C":
            return 0, self.imgSize
        if char == "G":
            return self.imgSize, self.imgSize
        if char == "T":
            return self.imgSize, 0
        else:
            return 0, 0
    
    def nextPoint(self, prevPoint, char):
        char_point = self.charPoint(char)
        return [int((prevPoint[0] + char_point[0])/2.0),
                int((prevPoint[1] + char_point[1])/2.0)]

    def getPoints(self, sequence):
        points = [[0, 0]]
        for s in sequence:
            points.append(self.nextPoint(points[(len(points)-1)], s))
        return points

    def pointsToGrid(self, points):
        grid = np.zeros([self.imgSize, self.imgSize], dtype=np.int)
        for point in points:
            grid[point[0], point[1]] = 255
        return grid

    def plot(self, points):
        #df = pd.DataFrame(points)
        #plt.scatter([p[0] for p in points], [p[1] for p in points])

        img = Image.new('1', (self.imgSize, self.imgSize), 'white')
        pixels = img.load()  # create pixel map
        grid = self.pointsToGrid(points)

    #    for i in range(img.size[0]):
    #        for j in range(img.size[1]):
    #            #pixels[i,j] = random.choice([0, 1])
    #            pixels[i,j] = grid[i, j]
        plt.axis('off')
        plt.imshow(grid, cmap='Greys', interpolation='nearest')
        plt.savefig("tmp/cgr_%d.png" % self.seqID, bbox_inches='tight')

    def run(self, seq, idx):
        sequence = list(str(seq.seq))
        output = open("tmp/atn_%d.pkl" % idx, 'wb')
        pickle.dump(seq.annotations, output)
        output.close()
        self.seqID = idx
        self.points = self.getPoints(sequence)
        self.plot(self.points)

    def rawpoints(self):
        return self.points


if __name__ == '__main__':
    #cgr = CGR()
    #cgr.run(list("CATTGTTGAGATCACATAATAATTGATCGAGTTAATCTGGAGGATCTGTTTACTTTGGTCACCCATGGGCATTTGCTGTTGAAGTGACCTAGATTTGCCATCGAGCCTCCTTGGGAGCTTTCTTGTTGGCGAGATCTAAACCCCTGCCCGGCGGAGTTGGGCGCCAAGTCATATGACACATAATTGGTGAAGGGGGTGGTAATCCTGCCCTGACCCTCCCCAAATTATTTTTTTAACAACTCTCAGCAACGGATATCTCGGCTCTTGCATCGATGAAGAACGCAGCGAAATGCGATAATGGTGTGAATTGCAGAATCCCGTGAACATCGAGTCTTTGAACGCAAGTTGCGCCCGAGGCCATCAGGCCAAGGGCACGCCTGCCTGGGCATTGCGAGTCATATCTCTCCCTTAATGAGGCTGTCCATACATACTGTTCAGCCGGTGCGGATGTGAGTTTGGCCCCTTGTTCTTTGGTACGGGGGGTCTAAGAGCTGCATGGGCTTTGGATGGTCCTAAATACGGAAAGAGGTGGACGAACTATGCTACAACAAAATTGTTGTGCAAATGCCCCGGTTGGCCGTTTAGTTGGGCC"))
    seq_recs_med = SeqIO.parse("data/mito_short.gbff", "gb")
    cgr = CGR()
    idx = 0
    for seq in seq_recs_med:
        idx += 1
        cgr.run(seq, idx)
