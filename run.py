from Order import Genomic


def main():
    path = "GRCh38_latest_genomic.fna"
    genomic = Genomic(path)

    ret = genomic.get_gene(6, 100000, 199999)
    order_0, _ = genomic.get_order(order=0)
    order_1, _ = genomic.get_order(order=1)
    order_2, _ = genomic.get_order(order=2)
    print({
        "order 0": order_0,
        "order 1": order_1,
        "order 2": order_2
    })

    prob, path = genomic.get_hmm_prob()

    print("prob: ", prob)

if __name__ == '__main__':
    main()