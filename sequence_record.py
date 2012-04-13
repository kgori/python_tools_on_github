#/usr/bin/env python
class Unaligned_Sequence_Record(object):
    """ Class for holding a unaligned sequences """
    def __init__(self, name, headers, sequences):
        self.name = name
        self.headers = headers
        self.sequences = sequences
        self.mapping = dict(zip(headers,sequences))
        self.get_alignment_columns()
        self.index = -1
        self.length = len(self.headers)
    def __iter__(self):
        return self
    def next(self):
        self.index += 1
        if self.index == self.length: 
            self.index = -1
            raise StopIteration
        return { 'header' : self.headers[self.index],
                 'sequence' : self.sequences[self.index]}
    def __str__(self):
        output_string = ""
        output_string += "Unaligned_Sequence_Record: {0}\n".format( self.name )
        for i in range(len(self.headers)):
            output_string += ">{0}\n{1}...({2})\n".format( self.headers[i], self.sequences[i][:50],len(self.sequences[i]) )
        output_string += "{0} sequences in record".format(len(self))
        return output_string
    def __len__(self):
        return self.length
    def get_alignment_columns(self):
        pass

    def sort_by_length(self):
        h, s = zip(*sorted(zip(self.headers,self.sequences), key = lambda item: len(item[1]),reverse=True))
        self.headers = h
        self.sequences = s

    def sort_by_name(self):
        h, s = zip(*sorted(zip(self.headers,self.sequences)))
        self.headers = h
        self.sequences = s

    def write_fasta(self,outfile="stdout",print_to_screen=False):
        lines = [">{0}\n{1}".format(h,seq) for h,seq in zip(self.headers,self.sequences)]
        s = '\n'.join(lines)
        s += '\n'
        if outfile == "stdout":
            print s
            return s
        elif outfile == "pipe":
            if print_to_screen: print s
            return s
        else: 
            open(outfile,'w').write(s)
            if print_to_screen: print s

    def write_nexus(self,outfile="stdout",sequence_type='protein'):
        maxlen = len(max(self.sequences,key=len))
        lines = [ "{0:<14} {1:-<{2}}".format(x,y,maxlen) for (x,y) in zip(self.headers,self.sequences) ]
        file_header =  "#NEXUS\n\n"
        file_header += "begin data;\n"
        file_header += "    dimensions ntax={0} nchar={1};\n".format(self.length,maxlen)
        file_header += "    format datatype={0} interleave=no gap=-;\n".format(sequence_type)
        file_header += "    matrix\n\n"

        file_footer = "    ;\n\nend;\n"

        s = file_header + '\n'.join(lines) + file_footer
        if outfile == "stdout":
            print s
            return s
        elif outfile == "pipe":
            return s
        else: open(outfile,'w').write(s)
    def write_phylip(self,outfile="stdout",sequence_type='protein',print_to_screen=False):
        maxlen = len(max(self.sequences,key=len))
        lines = [ "{0:<14} {1:-<{2}}".format(x,y,maxlen) for (x,y) in zip(self.headers,self.sequences) ]
        file_header = " {0} {1}\n".format(self.length,maxlen)
        s = file_header + '\n'.join(lines)
        s += '\n'
        if outfile == "stdout":
            print s
            return s
        elif outfile == "pipe":
            if print_to_screen: print s
            return s
        else: 
            open(outfile,'w').write(s)
            if print_to_screen: print s


class Aligned_Sequence_Record(Unaligned_Sequence_Record):
    """ Class for holding a sequence alignment """
    def __str__(self):
        output_string = ""
        output_string += "Aligned_Sequence_Record: {0}\n".format( self.name )
        for i in range(len(self.headers)):
            output_string += ">{0}\n{1}...({2})\n".format( self.headers[i], self.sequences[i][:50], len(self.sequences[i]) )
        output_string += "{0} sequences in record".format(len(self))
        return output_string
    def get_alignment_columns(self): 
        if not len(self.sequences) > 0:
            return
        self.columns = []
        for i in range(len(self.sequences[0])):
            column = ""
            for seq in self.sequences:
                column += seq[i]
            self.columns.append(column) 
    def sort_by_length(self):
        h, s = zip(*sorted(zip(self.headers,self.sequences), key = lambda item: len(item[1].replace('-','')),reverse=True))
        self.headers = h
        self.sequences = s
    
def get_fasta_file(fasta_file, name = "no name"):
    """ FASTA format parser: turns fasta file into Alignment_record object """
    headers = []
    sequences = []
    openfile = open(fasta_file,'r')

    #skip over file until first header is found
    while True:
        line = openfile.readline()
        if not line: return
        if line[0] == ">": break
        #we break the loop here at the start of the first record
    
    headers.append(line[1:].rstrip()) #chuck the first header into our list

    while True:
        line = openfile.readline()
        sequence_so_far = [] #build up sequence a line at a time in this list
        while True:
            if not line: break
            elif not line[0] == ">": 
                sequence_so_far.append(line.rstrip())
                line = openfile.readline()    
            else: break
        sequences.append("".join(sequence_so_far).replace(",",""))
        if not line: break
        headers.append(line[1:].rstrip())
    
    #check all sequences the same length
    first_seq_length = len(sequences[0])
    is_alignment = True
    for seq in sequences:
        if len(seq) != first_seq_length:
            is_alignment = False
            break
        else: continue

    #check same number of headers as sequences
    if len(headers) != len(sequences): print "Error matching all headers and sequences"

    if is_alignment: return Aligned_Sequence_Record(name, headers, sequences)
    else: return Unaligned_Sequence_Record(name, headers, sequences)

def get_phylip_file(phylip_file,name=None):
    """ PHYLIP format parser"""
    headers = []
    sequences = []
    openfile = open(phylip_file,'r')
    info = openfile.readline().split()
    num_taxa = info[0]
    seq_length = info[1]

    while True:
        line = openfile.readline().rstrip()
        if not line:
            break
        line = line.split()
        headers.append(line[0])
        sequences.append(line[1])

    return Aligned_Sequence_Record(name, headers, sequences)

def concatenate_2_alignments(alignment1, alignment2,name=None):
    try: assert set.intersection(set(alignment1.headers),set(alignment2.headers))==set(alignment1.headers)==set(alignment2.headers)
    except AssertionError: 
        print "Sequence labels do not match"
        return
    keys = alignment1.headers
    dict1 = dict(zip(alignment1.headers,alignment1.sequences))
    dict2 = dict(zip(alignment2.headers,alignment2.sequences))
    new_sequences = []
    for key in keys:
        new_sequences.append(dict1[key]+dict2[key])
    return Aligned_Sequence_Record(name,keys,new_sequences)

def concatenate_alignments(alignment_list,name=None):
    headers = alignment_list[0].headers
    data = [dict(zip(headers,alignment_list[0].sequences))]
    concatenation = ["" for x in headers]
    for i in range(1,len(alignment_list)):
        try: assert set.intersection(set(alignment_list[0].headers),set(alignment_list[i].headers))==set(alignment_list[0].headers)==set(alignment_list[i].headers)
        except AssertionError: 
            print "Sequence labels do not match between {0} and {1}".format(alignment_list[0],alignment_list[i])
            return
    
        data.append(dict(zip(alignment_list[i].headers,alignment_list[i].sequences)))
    
    for i in range(len(headers)):
        h = headers[i]
        for j in range(len(data)):
            concatenation[i]+=data[j][h]

    return Aligned_Sequence_Record(name,headers,concatenation)
