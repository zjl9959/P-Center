#!usr/bin/python3
import csv

# ALGORITHM_INDEX = 1
# INSTANCE_INDEX = 2
# VALID_INDEX = 7
# OBJ_INDEX = 8
# CONFIG_INDEX = 9

def _parser_log():
    res = {}
    count = 0
    with open('log.csv') as f:
        records = csv.reader(f)
        for record in records:  # search every log and record key information
            if record[7] == '1':
                cfg = record[1] + ';k=' + record[9][1:record[9].find('/')]
                instname = record[2][:record[2].find('.')] + '-' + record[3]
                if res.__contains__(cfg) is False:  # creat a new dictionary for new configuration
                    res[cfg] = {}
                if res[cfg].__contains__(instname):  # get average value for the same instance
                    res[cfg][instname] = (res[cfg][instname] + float(record[8])) * 0.5
                else:
                    res[cfg][instname] = float(record[8])
            count = count + 1
    print('Parsed %d records' %(count))
    return res


def _load_best_obj():
    res = {}
    with open('instance/best_result.csv') as f:
        infos = csv.reader(f)
        for info in infos:
            if info[0] != '-i':
                res[info[0] + '-' + info[2]] = float(info[3])
    print("Load %d best results." %(res.__len__()))
    return res


def main():
    result = _parser_log()
    best = _load_best_obj()
    for cfg in result:
        bigger = 0
        smaller = 0
        same = 0
        unknown = 0
        delt = 0.0
        print("%s has %d instance:" %(cfg, result[cfg].__len__()))
        for inst in result[cfg]:
            if best.__contains__(inst):
                if result[cfg][inst] > best[inst]:
                    bigger = bigger + 1
                elif result[cfg][inst] == best[inst]:
                    same = same + 1
                else:
                    smaller = smaller + 1
                delt = delt + (result[cfg][inst] - best[inst])
            else:
                unknown = unknown + 1
        print("<:%d\t=:%d\t>:%d\t?:%d\tdelt:%d"%(smaller, same, bigger, unknown, delt))


if __name__ == '__main__':
    main()