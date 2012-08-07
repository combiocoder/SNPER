# Classes to support snp mapping
class ExonInfo:
    """A simple class to store information about introns of a gene"""
    seq = ''
    def __init__(self,info_list):
        self.refseq_id = info_list[1]
        self.chromosome = info_list[2]
        self.strand = info_list[3]
        self.cds_start = int(info_list[6])
        self.cds_end = int(info_list[7])
        self.exon_count = int(info_list[8])
        self.exon_starts = map(int,info_list[9].split(',')[:-1])
        self.exon_ends = map(int,info_list[10].split(',')[:-1])
        self.exon_frames = map(int,info_list[15].split(',')[:-1])
        self.hit_exon_start = 1
        for i in range(self.exon_count):
            if self.cds_start >= self.exon_starts[i] and \
                                 self.cds_start < self.exon_ends[i]:
                self.cds_start_pos_in_exons = self.hit_exon_start + (
                                          self.cds_start - self.exon_starts[i])
            else: self.hit_exon_start += (self.exon_ends[i] - 
                                          self.exon_starts[i]) 
        self.hit_exon_start = 1
        for i in range(self.exon_count):
            if self.cds_end > self.exon_starts[i] and \
                                 self.cds_end <= self.exon_ends[i]:
                self.cds_end_pos_in_exons = self.hit_exon_start + (
                                          self.cds_end - self.exon_starts[i])
            else: self.hit_exon_start += (self.exon_ends[i] - 
                                          self.exon_starts[i]) 
