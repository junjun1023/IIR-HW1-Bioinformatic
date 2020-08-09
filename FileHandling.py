

path = "GRCh38_latest_genomic.fna"

def get_gene(path, ith, start=0, end=-1):
    flag = False
    i = 0

    ret = ""
    ret_cursor = 0

    cursor = 0
    r_flag = False

    with open(path) as f_obj:
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
 
    return ret


ret = get_gene(path, 6, start=100000, end=199999)

print(ret)

    

