import itertools
import math

class Genomic(object):

    def __init__(self, path):
        super(Genomic, self).__init__()
        self._path = path


    def get_gene(self, ith, start=0, end=-1):
        flag = False
        i = 0

        ret = ""
        ret_cursor = 0

        cursor = 0
        r_flag = False

        with open(self._path) as f_obj:
            for line in f_obj:

                if "Homo sapiens chromosome {ith}".format(ith=ith) in line:
                    flag = True
                    continue
                elif "Homo sapiens chromosome {ith}".format(ith=ith+1) in line:
                    break
                elif flag == False:
                    continue

                line = line.replace("\n", "")
                length = len(line)

                if cursor == start:
                    r_flag = True
                elif end == -1:
                    r_flag = True
                elif cursor > end:
                    break

            
                if r_flag:
                    if end == -1:
                        ret = ret + line
                    if cursor+length > end:
                        need = end - cursor + 1
                        ret = ret + line[:need]
                        break  
                    else:
                        ret = ret + line

                cursor = cursor + length
 
        self._source = ret

        return ret



    def get_order(self, source, order=0, base='atcg'):

        key_value = {}
        keys = itertools.product(base, repeat=order+1)
        for key in keys:
            k = ''.join(key)
            key_value[k] = 0

        if source == None:
            if self._source == None:
                print('Error:- Missing source, execute function -- get_gene first')
            else:
                source = self._source

        source = source.lower()
        total = len(source) - order
        for i in range(len(source)-order):
            k = source[i:i+order+1]
            key_value[k] += 1

        prob_key_value = {}
        for k, v in key_value.items():
            prob = v/total
            prob_ln = math.log2(prob)
            prob_key_value[k] = prob_ln

        prob = 0
        for k, v in prob_key_value.items():
            times = key_value[k]
            prob = prob + times * v

        # print(prob)
            

        return prob
        

    