import itertools
import math

class Genomic(object):

    _start = {
                's1': 0.3,
                's2': 0.3,
                's3': 0.4
            }
    _transmission = {
                's1': {
                    's1': 0.5,
                    's2': 0.25,
                    's3': 0.25
                },
                's2': {
                    's1': 0.25,
                    's2': 0.5,
                    's3': 0.25
                },
                's3': {
                    's1': 0.25,
                    's2': 0.25,
                    's3': 0.5
                }
            }
    _emission = {
                's1': {
                    'a': 0.4,
                    't': 0.2,
                    'c': 0.2,
                    'g': 0.2
                },
                's2': {
                    'a': 0.2,
                    't': 0.4,
                    'c': 0.2,
                    'g': 0.2
                },
                's3': {
                    'a': 0.2,
                    't': 0.2,
                    'c': 0.4,
                    'g': 0.2
                }
            }
    _states = ('s1', 's2', 's3')


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



    def get_order(self, order=0, base='atcg', source=None):

        key_value = {}

        for i in range(order+1):
            keys = itertools.product(base, repeat=i+1)
            for key in keys:
                k = ''.join(key)
                key_value[k] = 0

        if source == None:
            if self._source == None:
                print('Error:- Missing source, execute function -- get_gene first')
                return
            else:
                source = self._source

        source = source.lower()
        total = len(source)

        for i in range(total):
            for j in range(order+1):
                k = source[i:i+j+1]
                if k in key_value:
                    key_value[k] += 1

        # print(key_value)

        prob_key_value = {}
        for k, v in key_value.items():
            prob = 2
            if len(k) == 1:
                prob = v/total
            else:
                k2 = k[:-1]
                
                if k2 in key_value:
                    prob = v/key_value[k2]
                    # print(k, k2, prob)
                    
            prob_ln = math.log2(prob)
            prob_key_value[k] = prob_ln

        # print(prob_key_value)
        
        prob = 0
        for i in range(order):
            key = source[0:i+1]
            prob += prob_key_value[key]

        for i in range(total-order):
            key = source[i:i+order+1]
            prob += prob_key_value[key]


        # print(prob)
            

        return prob, key_value
        

    
    def get_hmm_prob(self, states=None, start=None, transmission=None, emission=None, source=None):
        # 在這個狀態往下一個狀態看
        if source == None:
            if self._source == None:
                print('Error:- Missing source, execute function -- get_gene first')
                return
            else:
                source = self._source
        source = source.lower()

        if states==None:
            states = self._states

        if start==None:
            start = self._start

        if transmission==None:
            transmission = self._transmission
        
        if emission==None:
            emission = self._emission

        path = {}
        V = [{}]


        # Initialize base cases (t == 0)
        for y in states:
            V[0][y] = math.log2(start[y]) + math.log2(emission[y][source[0]])
            path[y] = [y]


        # Run Viterbi for t > 0
        for t in range(1,len(source)):
            V.append({})
            newpath = {}

            for y in states:

                (prob, state) = max( [(V[t-1][y0] + math.log2(transmission[y0][y]) + math.log2(emission[y][source[t]]), y0)  for y0 in states])
                # print(prob, state)
                
                V[t][y] = prob

                # print(path)
                # print(state)
                # print(path[state])
                newpath[y] = path[state] + [y]

                # Don't need to remember the old paths
            path = newpath

        (prob, state) = max([(V[len(source) - 1][y], y) for y in states])

        fp = open("state_path.txt", "w")
        fp.write(str(path[state]))
        fp.close()

        return (prob, path[state])









    # def _source_valid(self, source):
    #     if source == None:
    #         if self._source == None:
    #             self.
    #             return 
    #         else:
    #             source = self._source
    #             return source
    #     return 
    