import os
import psutil
import time
from Order import Genomic

current_process = psutil.Process(os.getpid())
f = open("memory_usage.txt", "w")

def main():
    path = "GRCh38_latest_genomic.fna"
    genomic = Genomic(path)

    start_time = time.time()
    start_mem = current_process.memory_info().rss # in byte
    ret = genomic.get_gene(6, 100000, 199999)
    end_mem = current_process.memory_info().rss
    end_time = time.time()
    f.write("Done Extract: " + str(end_mem-start_mem) + " bytes \n")
    f.write("Done Extract: " + str(end_time-start_time) + " sec \n")

    start_time = time.time()
    start_mem = current_process.memory_info().rss # in byte
    order_0, _ = genomic.get_order(order=0)
    end_mem = current_process.memory_info().rss
    end_time = time.time()
    f.write("Done Order 0: " + str(end_mem-start_mem) + " bytes \n")
    f.write("Done Order 0: " + str(end_time-start_time) + " sec \n")

    start_time = time.time()
    start_mem = current_process.memory_info().rss # in byte
    order_1, _ = genomic.get_order(order=1)
    end_mem = current_process.memory_info().rss
    end_time = time.time()
    f.write("Done Order 1: " + str(end_mem-start_mem) + " bytes \n")
    f.write("Done Order 1: " + str(end_time-start_time) + " sec \n")

    start_time = time.time()
    start_mem = current_process.memory_info().rss # in byte
    order_2, _ = genomic.get_order(order=2)
    end_mem = current_process.memory_info().rss
    end_time = time.time()
    f.write("Done Order 2: " + str(end_mem-start_mem) + " bytes \n")
    f.write("Done Order 2: " + str(end_time-start_time) + " sec \n")

    print({
        "order 0": order_0,
        "order 1": order_1,
        "order 2": order_2
    })

    start_time = time.time()
    start_mem = current_process.memory_info().rss # in byte
    prob, path = genomic.get_hmm_prob()
    end_mem = current_process.memory_info().rss
    end_time = time.time()
    f.write("Done HMM: " + str(end_mem-start_mem) + " bytes \n")
    f.write("Done HMM: " + str(end_time-start_time) + " sec \n")

    print("prob: ", prob)

if __name__ == '__main__':

    main()

    f.close()