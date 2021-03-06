#!usr/bin/python3
import os
import csv

configs = list()

def _load_configs():
    global configs
    with open(os.getcwd() + '/config.txt') as f:
        labels = f.readline()[:-1].split(' ')
        lines = f.readlines()  # parser infomation
        for line in lines:
            if line == 'break\n':
                break
            info = line[:-1].split(' ')
            cfg = list()
            for i in range(len(labels)):
                cfg.append((labels[i],info[i]))
            configs.append(cfg)


def _benchmark_test():
    global configs
    for cfg in configs:
        cmd = os.getcwd() + '/PCenter.exe '
        for switch in cfg:  # creat shell commend
            if switch[0] == '-i' or switch[0] == '-p':
                cmd += switch[0] + ' '
                if switch[0] == '-i' and switch[1].count('pmed'):
                    cmd += switch[1] + '.txt' + ' '
                elif switch[0] == '-i':
                    cmd += switch[1] + '.tsp' + ' '
                else:
                    cmd += switch[1] + ' '
        os.system(cmd)
        

def main():
    _load_configs()
    _benchmark_test()


if __name__ == '__main__':
    main()
