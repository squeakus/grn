# This file is part of Python GRN implementation.
# Architype is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# Architype is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Lesser General Public License for more details.
# You should have received a copy of the GNU Lesser General Public License 
# along with GRN.  If not, see <http://www.gnu.org/licenses/>.
# Author Jonathan Byrne 2014

"""
This class moves the promoter site to the top of the gene location
"""

def reorder_gene(index, seq, genome):
    """chop gene into segments and reorder:
    promoter + enhance + inhibit + gene_info"""
    offset = 96 - len(seq)
    start_index = index - offset
    gene = genome[start_index:start_index+256]
    enh = gene[0:32]
    inh = gene[32:64]
    prom = gene[64:96]
    info = gene[96:256]
    flipped_promoter = prom[32-len(seq):32] + prom[0:32-len(seq)]
    new_gene = flipped_promoter + enh + inh + info

    #check against miguels code
    # enh_str = [str(x) for x in enh]
    # enh_str = ''.join(enh_str)
    # inh_str = [str(x) for x in inh]
    # inh_str = ''.join(inh_str)
    # prom_str = [str(x) for x in prom]
    # prom_str = ''.join(prom_str)
    # print "en", int(enh_str,2), "in", int(inh_str,2), "pro", int(prom_str,2)
    return new_gene

def main():
    """check it is finding gene indexes correctly"""
    sequences = [[0] * 8, [1] * 8]
    index = 0
    gen_file = open("miguel_parsed.txt",'r')
    genome = eval(gen_file.readline())
    while index < len(genome) -256:
        for seq in sequences:
            #print index, "comparing:", genome[index:index+len(seq)], "and", seq
            found = genome[index:index+len(seq)] == seq
            if found:
                print "found!", genome[index:index+len(seq)], "at", index
                reorder_gene(index, seq, genome)
                index += 255
        index += 1

if __name__ == "__main__":
    main()
