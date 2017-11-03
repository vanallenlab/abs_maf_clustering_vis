from abs_maf_analyzer.abs_maf_analyzer import AbsMafAnalyzer


def main():
    path = '/Users/erickofman/Documents/Projects/OralCancer/Copy Number Alteration Data/OCSCC-OC011-TP-NB-SM-F3R74-SM-F3R84_ABS_MAF.txt'
    #path = '/Users/erickofman/Documents/Projects/OralCancer/Copy Number Alteration Data/OCSCC-OC020-TP-NB-SM-F3R7A-SM-F3R88_ABS_MAF.txt'
    #path = '/Users/erickofman/Documents/Projects/OralCancer/Copy Number Alteration Data/OCSCC-OC026-TP-NB-SM-F3R7G-SM-F3R8A_ABS_MAF.txt'

    a = AbsMafAnalyzer(path, accession="OC")
    a.cluster()

if __name__ == '__main__':
    main()