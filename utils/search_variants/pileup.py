# -*- coding: utf8 -*-

import pileup_unit
import settings
import common
from common import log, CustomError


class PileupLine(object):
    def __init__(self, chromosome, position, ref, depth, bases, qualities):
        self.chromosome = chromosome
        self.position = position
        self.ref = ref
        self.depth = depth
        self.qualities = qualities
        self.bases = bases

    def __init__(self, line):
        if line == u'':
            raise IOError
        r = line.split(u"\t")
        self.chromosome = r[0]
        self.position = int(r[1])
        self.ref = r[2]
        self.depth = int(r[3])
        self.bases = r[4]
        self.qualities = r[5].rstrip()

    def get_bases_with_qualities(self):
        data = []
        # remove the beginning mark
        i = 0 # for bases
        j = 0 # for qualities
        while i < len(self.bases):
            # log.debug("{0}, {1}: {2}, {3}".format(i, j, self.bases[i], self.qualities[j]))
            if self.bases[i] == u".":
                data.append(pileup_unit.Base(self.ref, 1, self.qualities[j]))
                i += 1
                j += 1
            elif self.bases[i] == u",":
                data.append(pileup_unit.Base(self.ref, 0, self.qualities[j]))
                i += 1
                j += 1
            elif self.bases[i] == u"$":
                i += 1
            elif self.bases[i] == u"^":
                i += 2  # with map quality
                continue
            elif self.bases[i] in u"+-":
                # read the length
                tmp_i = i + 1
                length_str = u""
                while True:
                    try:
                        int(self.bases[tmp_i])
                        length_str += self.bases[tmp_i]
                        tmp_i += 1
                    except ValueError:
                        break
                #log.debug("length_str: {0}".format(length_str))
                seq_length = int(length_str)
                # get the inserted or deleted sequence
                seq = self.bases[tmp_i:tmp_i + seq_length]
                if seq[0] in u"ATGCN":
                    strand = 1
                else:
                    strand = 0
                if self.bases[i] == u"+":
                    data.append(pileup_unit.Ins(seq, strand, None))
                else:
                    data.append(pileup_unit.Del(seq, strand, None))
                i = tmp_i + seq_length
            elif self.bases[i] in u"ATGCN":
                data.append(pileup_unit.Base(self.bases[i], 1, self.qualities[j]))
                i += 1
                j += 1
            elif self.bases[i] in u"atgcn":
                data.append(pileup_unit.Base(self.bases[i].upper(), 0, self.qualities[j]))
                i += 1
                j += 1
            elif self.bases[i] == u"*":
                # deletion
                i += 1
                j += 1
            else:
                log.warning(u"i={0}, self={1}".format(i, self))
                raise Exception(u"something_wrong")

        if i != len(self.bases) or j != len(self.qualities):
            log.error(u"i:{0}={1} , j:{2}={3}".format(i, len(self.bases), j, len(self.qualities)))
            log.error(self)
            raise Exception(u"did_not_reach_the_lase_base")
        return data

    def __str__(self):
        return u"{0}:{1},{2},depth=({3}),bases=({4}),qualities=({5})".format(self.chromosome, self.position, self.ref, self.depth,
                                                                             self.bases, self.qualities)

    def __repr__(self):
        return self.__str__()
