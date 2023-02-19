import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
import os
import json

class Histogram:
    def __init__(self, data_file, random_file, cfg_file, angular_boxes):
        self.data_file = data_file[:-5]
        self.rand_file = random_file[:-5]
        self.cfg_file = cfg_file
        self.partitions = angular_boxes

        with open('./params/'+cfg_file) as cfgDict:
            cfg_set = json.load(cfgDict)
            self.sMin = cfg_set['sMin']
            self.sMax = cfg_set['sMax']
            self.sBins = cfg_set['sBinN']
            

    
    """

    @param: number of partitions in the catalog
    @return: w0 and b0 in that order
    
    """
    def process_plain_grid(self, partitions):
        w0 = []
        b0 = []

        hdul_data = fits.open("./data/{i}.fits".format(i=self.data_file))
        hdul_randoms = fits.open("./data/{i}.fits".format(i=self.rand_file))
        data = hdul_data[1].data.field(3)
        randoms = hdul_randoms[1].data.field(3)

        nD = np.sum(data)
        nR = np.sum(randoms)
        
        for i in range(1,self.partitions+1):
            hdul = fits.open("./grids/{h}_ex_par/W_p{i}_of_{m}.fits".format(h=self.data_file, i=i, m=self.partitions))
            data = hdul[1].data.field(0)
            w0 = np.concatenate((w0, data), axis=None) 

        for i in range(1,self.partitions+1):
            hdul = fits.open("./grids/{h}_ex_par/B_p{i}_of_{m}.fits".format(h=self.rand_file,i=i, m=self.partitions))
            data = hdul[1].data.field(0)
            b0 = np.concatenate((b0, data), axis=None)

        b0 = b0*(nD/nR)
        return w0, b0 
        

    """
    @param: number of partitions, number of sbins
    @return: dictionary of concatenated W1 values
    """
    def process_all_w(self, partitions, sbins):
        w1Bins = {}  # dict that stores info for all s bins
        for s in range(sbins) :
            this_sbin = "w1_s{s}".format(s=s)
            oneBin = []

            for i in range(1,partitions+1):
                hdul = fits.open("./grids/{h}_ex_par/W_{k}_0_0_RE_p{i}_of_{m}.fits".format(h=self.data_file,k=s,i=i, m=partitions))
                data = hdul[1].data.field(0)
                oneBin = np.concatenate((oneBin, data), axis=None)
            
            w1Bins[this_sbin] = oneBin

        return w1Bins



    """
    @param: number of partitions, number of sbins
    @return: dictionary of concatenated B1 values

    """
    def process_all_b(self, partitions, sbins):

        # calculate normalizing coefficient
        hdul_data = fits.open("./data/{i}.fits".format(i=self.data_file))
        hdul_randoms = fits.open("./data/{i}.fits".format(i=self.rand_file))
        data = hdul_data[1].data.field(3)
        randoms = hdul_randoms[1].data.field(3)

        nD = np.sum(data)
        nR = np.sum(randoms)  # do this for B0 as well

        b1Bins = {}  # dict that stores info for all s bins
        for s in range(sbins) :
            this_sbin = "b1_s{s}".format(s=s)
            oneBin = []

            for i in range(1,partitions+1):
                hdul = fits.open("./grids/{h}_ex_par/B_{k}_0_0_RE_p{i}_of_{m}.fits".format(h=self.rand_file,k=s,i=i,m=partitions))
                data = hdul[1].data.field(0)
                oneBin = np.concatenate((oneBin, data), axis=None)

            b1Bins[this_sbin] = oneBin*(nD/nR)  

        return b1Bins


    """

    @param: s_index to plot, w0 from process_plain_grid, b0 from process_plain_grid, w1_dict from process_all_w
    @return: none

    """
    def plot_hist(self, s_index, w0, b0, w1_dict, b1_dict):
        toPlot_w1 = "w1_s{i}".format(i=s_index)
        toPlot_b1 = "b1_s{j}".format(j=s_index)
        w1 = w1_dict[toPlot_w1]
        b1 = b1_dict[toPlot_b1]

        w0w1 = w0 * w1
        b0b1 = b0 * b1
        b0b1_mean = np.mean(b0b1)

        normalized = w0w1/b0b1_mean

        if not os.path.exists("./hist_plots"):    
            os.makedirs("./hist_plots")

        if not os.path.exists("./hist"):    
            os.makedirs("./hist")

        s_MPCH = s_index * (self.sMax - self.sMin)/(self.sBins-1) + self.sMin
        lbl = "s = " + str(s_MPCH)
        plt.figure()
        n, bins, patches = plt.hist(normalized, 100, histtype='stepfilled', alpha=0.25, color='greenyellow', edgecolor='black', lw=2, label=lbl + r" $h^{-1}$Mpc")
        # plt.hist(normalized, 100, histtype='stepfilled', edgecolor='black', facecolor="None")
        plt.title("{t}".format(t=self.data_file))
        plt.yscale("log")
        plt.xlabel("W0W1/<B0B1>")
        plt.ylabel("pdf")        
        plt.legend()
        plt.savefig("./hist_plots/{f}_s{s}.png".format(f=self.data_file, s=s_index))

        f = open("./hist/{f}_s{s}.txt".format(f=self.data_file, s=s_index), "w")
        f.write("count\n")
        for x in n:
            f.write(str(x) + "\n")
        f.write("\n")
        
        f.write("bins\n")
        for y in bins:
            f.write(str(y) + "\n")
        
        f.close()

        return n, bins, patches



def ProcessHistogram(data_file:str, random_file:str, config_file: str, angular_boxes: int):
    histprocessor = Histogram(data_file, random_file, config_file, 45)

    print("Processing plain grids...\n")
    w0, b0 = histprocessor.process_plain_grid(histprocessor.partitions)
    print("Plain grids processed\n")
    print("Processing data catalog...\n")
    w1_dict = histprocessor.process_all_w(histprocessor.partitions,histprocessor.sBins)
    print("Data catalog processed\n")
    print("Processing random catalog...\n")
    b1_dict = histprocessor.process_all_b(histprocessor.partitions,histprocessor.sBins)
    print("Random catalog processed\n")
    print("Saving histogram figures...\n")

    for s in range(histprocessor.sBins):
        histprocessor.plot_hist(s, w0, b0, w1_dict, b1_dict)
    print("Done!\n")


if __name__ == "__main__":
    ProcessHistogram('fnl0sim1.fits', 'randoms.fits', 'ex_par.txt', 45)
