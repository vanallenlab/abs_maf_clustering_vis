The AbsMafAnalysis workflow can be used to visualize the data from an *ABS_MAF.txt absolute file output.
Both detection power and genomic locus are plotted for all SNPs, along with clustering information based
on cancer cell fraction.

Example usage in a script:

    from abs_maf_analyzer.abs_maf_analyzer import AbsMafAnalyzer

    def main():
        path = '/Users/erickofman/Documents/Projects/OralCancer/OCSCC-OC011-TP-NB-SM-F3R74-SM-F3R84_ABS_MAF.txt'
        a = AbsMafAnalyzer(path, accession="OC011")

        # k is not specified -- all values of k between 1 and 10 will be tested for optimal number of clusters. This
        # takes around 12 seconds.
        a.cluster()

        path = '/Users/erickofman/Documents/Projects/OralCancer/OCSCC-OC011-TP-NB-SM-F3R74-SM-F3R84_ABS_MAF.txt'
        b = AbsMafAnalyzer(path, accession="OC011")

        # k is specified -- value specified will be used for number of clusters
        b.cluster(k=5)

    if __name__ == '__main__':
        main()

* The k-means implementation used here, which leverages pandas, is loosely based on an implementation by Jack Maney,
which can be found here: https://github.com/jackmaney/k-means-plus-plus-pandas