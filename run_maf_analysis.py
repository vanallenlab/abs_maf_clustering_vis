from abs_maf_analyzer.abs_maf_analyzer import AbsMafAnalyzer


def main():
    path = '/Users/erickofman/Documents/Projects/OralCancer/Copy Number Alteration Data/OCSCC-OC011-TP-NB-SM-F3R74-SM-F3R84_ABS_MAF.txt'
    a = AbsMafAnalyzer(path, accession="OC011")
    a.cluster()

    path = '/Users/erickofman/Documents/Projects/OralCancer/Copy Number Alteration Data/OCSCC-OC020-TP-NB-SM-F3R7A-SM-F3R88_ABS_MAF.txt'
    a = AbsMafAnalyzer(path, accession="OC020")
    a.cluster()

    path = '/Users/erickofman/Documents/Projects/OralCancer/Copy Number Alteration Data/OCSCC-OC026-TP-NB-SM-F3R7G-SM-F3R8A_ABS_MAF.txt'
    a = AbsMafAnalyzer(path, accession="OC026")
    a.cluster()

    path = '/Users/erickofman/Documents/Projects/OralCancer/Copy Number Alteration Data/OCSCC-OC033-TP-NB-SM-F3R7J-SM-F3R8D_ABS_MAF.txt'
    a = AbsMafAnalyzer(path, accession="OC033")
    a.cluster()

    path = '/Users/erickofman/Documents/Projects/OralCancer/Copy Number Alteration Data/OCSCC-OC009-TP-NB-SM-F3R72-SM-F3R82_ABS_MAF.txt'
    a = AbsMafAnalyzer(path, accession="OC009")
    a.cluster()

    path = '/Users/erickofman/Documents/Projects/OralCancer/Copy Number Alteration Data/OCSCC-OC031-TP-NB-SM-F3R7I-SM-F3R8C_ABS_MAF.txt'
    a = AbsMafAnalyzer(path, accession="OC031")
    a.cluster()


if __name__ == '__main__':
    main()
