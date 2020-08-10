from Order import Genomic


def main():
    path = "GRCh38_latest_genomic.fna"
    genomic = Genomic(path)

    ret = genomic.get_gene(6, 100000, 199999)
    order_0 = genomic.get_order(None, order=0)
    order_1 = genomic.get_order(None, order=1)
    order_2 = genomic.get_order(None, order=2)
    print({
        "order 0": order_0,
        "order 1": order_1,
        "order 2": order_2
    })

if __name__ == '__main__':
    main()