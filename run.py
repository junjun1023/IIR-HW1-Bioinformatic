import os
import psutil
from Order import Genomic

current_process = psutil.Process(os.getpid())
f = open("memory_usage.txt", "a")

def main():
    path = "GRCh38_latest_genomic.fna"
    genomic = Genomic(path)

    ret = genomic.get_gene(6, 100000, 199999)
    end_mem = current_process.memory_info().rss
    f.write("Done Extract: " + str(end_mem) + "\n")

    order_0, _ = genomic.get_order(order=0)
    end_mem = current_process.memory_info().rss
    f.write("Done Order 0: " + str(end_mem) + "\n")

    order_1, _ = genomic.get_order(order=1)
    end_mem = current_process.memory_info().rss
    f.write("Done Order 1: " + str(end_mem) + "\n")

    order_2, _ = genomic.get_order(order=2)
    end_mem = current_process.memory_info().rss
    f.write("Done Order 2: " + str(end_mem) + "\n")

    print({
        "order 0": order_0,
        "order 1": order_1,
        "order 2": order_2
    })

    prob, path = genomic.get_hmm_prob()
    end_mem = current_process.memory_info().rss
    f.write("Done HMM: " + str(end_mem) + "\n")

    print("prob: ", prob)

if __name__ == '__main__':

    start_mem = current_process.memory_info().rss # in byte
    f.write("main:- start " + str(start_mem) + "\n")

    main()

    f.close()