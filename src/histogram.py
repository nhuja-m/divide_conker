import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from scipy.stats import moment
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
    def process_all_w(self, partitions):
        w1Bins = {}  # dict that stores info for all s bins
        for s in range(self.sBins) :
            this_sbin = s
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
    def process_all_b(self, partitions):

        # calculate normalizing coefficient
        hdul_data = fits.open("./data/{i}.fits".format(i=self.data_file))
        hdul_randoms = fits.open("./data/{i}.fits".format(i=self.rand_file))
        data = hdul_data[1].data.field(3)
        randoms = hdul_randoms[1].data.field(3)

        nD = np.sum(data)
        nR = np.sum(randoms)  # do this for B0 as well

        b1Bins = {}  # dict that stores info for all s bins
        for s in range(self.sBins) :
            this_sbin = s
            oneBin = []

            for i in range(1,partitions+1):
                hdul = fits.open("./grids/{h}_ex_par/B_{k}_0_0_RE_p{i}_of_{m}.fits".format(h=self.rand_file,k=s,i=i,m=partitions))
                data = hdul[1].data.field(0)
                oneBin = np.concatenate((oneBin, data), axis=None)

            b1Bins[this_sbin] = oneBin*(nD/nR)  

        return b1Bins
    
    """
    Combines all the information to calculate required parameter and stores internally in master_info dictionary
    """
    
    def combine(self, w0, b0, w1_dict, b1_dict):
        combined_dict = {} # dict that combines info for all s bins
        for s in range(self.sBins):
            w1 = w1_dict[s]
            b1 = b1_dict[s]

            w0w1 = w0 * w1
            b0b1 = b0 * b1
            b0b1_mean = np.mean(b0b1)

            combined_dict[s] = w0w1/b0b1_mean
        
        self.master_info = combined_dict


    # Useful to calculate up to 7 moments
    # Stores output in folder called hist
    def calc_moment(self):

        if not os.path.exists("./hist"):    
            os.makedirs("./hist")

        moment1 = [None] * self.sBins
        moment2 = [None] * self.sBins
        moment3 = [None] * self.sBins
        moment4 = [None] * self.sBins
        moment5 = [None] * self.sBins
        moment6 = [None] * self.sBins
        moment7 = [None] * self.sBins


        for s in range(self.sBins):
            normalized = self.master_info[s]
            moment1[s] = np.mean(normalized)
            moment2[s] = moment(normalized, moment=2)
            moment3[s] = moment(normalized, moment=3)
            moment4[s] = moment(normalized, moment=4)
            moment5[s] = moment(normalized, moment=5)
            moment6[s] = moment(normalized, moment=6)
            moment7[s] = moment(normalized, moment=7)


        self.m1 = moment1
        self.m2 = moment2
        self.m3 = moment3
        self.m4 = moment4
        self.m5 = moment5
        self.m6 = moment6
        self.m7 = moment7


        f = open("./hist/{f}_moments.txt".format(f=self.data_file), "w")
        
        f.write(self.data_file + "\n")

        m1 = ""
        for x in moment1:
            m1 += str(x) + " "
        f.write(m1 + "\n")

        m2 = ""
        for x in moment2:
            m2 += str(x) + " "
        f.write(m2 + "\n")

        m3 = ""
        for x in moment3:
            m3 += str(x) + " "
        f.write(m3 + "\n")

        m4 = ""
        for x in moment4:
            m4 += str(x) + " "
        f.write(m4 + "\n")

        m5 = ""
        for x in moment5:
            m5 += str(x) + " "
        f.write(m5 + "\n")
        
        m6 = ""
        for x in moment6:
            m6 += str(x) + " "
        f.write(m6 + "\n")
        
        m7 = ""
        for x in moment7:
            m7 += str(x) + " "
        f.write(m7 + "\n")
        
        f.close()


def ProcessHistogram(data_file:str, random_file:str, config_file: str, angular_boxes: int):
    histprocessor = Histogram(data_file, random_file, config_file, 45)
    print("Working on {d}".format(d=data_file))
    print("Processing plain grids...\n")
    w0, b0 = histprocessor.process_plain_grid(histprocessor.partitions)
    print("Plain grids processed\n")
    print("Processing data catalog...\n")
    w1_dict = histprocessor.process_all_w(histprocessor.partitions)
    print("Data catalog processed\n")
    print("Processing random catalog...\n")
    b1_dict = histprocessor.process_all_b(histprocessor.partitions)
    print("Random catalog processed\n")
    print("Combining calculated data...\n")
    histprocessor.combine(w0, b0, w1_dict, b1_dict)
    print("Done!\n")

    return histprocessor

# def PlotHistograms(fnl0Object:Histogram, fnl100Object:Histogram):
#     if not os.path.exists("./hist_plots"):    
#             os.makedirs("./hist_plots")

#     if not os.path.exists("./hist"):    
#         os.makedirs("./hist")

#     if not os.path.exists("./hist/fnl0"):    
#         os.makedirs("./hist/fnl0")

#     if not os.path.exists("./hist/fnl100"):    
#         os.makedirs("./hist/fnl100")

#     for s_index in range(fnl0Object.sBins):
#         s_MPCH = s_index * (fnl0Object.sMax - fnl0Object.sMin)/(fnl0Object.sBins-1) + fnl0Object.sMin
#         lbl = "s = " + str(s_MPCH)

#         plt.figure()
#         n0, bins0, patches0 = plt.hist(fnl0Object.master_info[s_index], 100, histtype='stepfilled', alpha=0.50, color='greenyellow', edgecolor='black', label="fnl 0")
#         n100, bins100, patches100 = plt.hist(fnl100Object.master_info[s_index], 100, histtype='stepfilled', alpha=0.50, color='skyblue', edgecolor='black', label="fnl 100") # plt.hist(normalized, 100, histtype='stepfilled', edgecolor='black', facecolor="None")
        
#         plt.title("{t} and {u} " +  lbl + r" $h^{-1}$Mpc".format(t=fnl0Object.data_file, u = fnl100Object.data_file))
#         plt.yscale("log")
#         plt.xlabel("W0W1/<B0B1>")
#         plt.ylabel("count log scale")        
#         plt.legend()
#         plt.savefig("./hist_plots/sim{x}_s{s}.png".format(x=fnl0Object.data_file[7:], s=s_index))

#         f = open("./hist/fnl0/{f}_s{s}.txt".format(f=fnl0Object.data_file, s=s_index), "w")
#         f.write("count\n")
#         for x in n0:
#             f.write(str(x) + "\n")
#         f.write("\n")
        
#         f.write("bins\n")
#         for y in bins0:
#             f.write(str(y) + "\n")
        
#         f.close()

#         g = open("./hist/fnl100/{f}_s{s}.txt".format(f=fnl100Object.data_file, s=s_index), "w")
#         g.write("count\n")
#         for x in n100:
#             g.write(str(x) + "\n")
#         g.write("\n")
        
#         g.write("bins\n")
#         for y in bins100:
#             g.write(str(y) + "\n")
        
#         g.close()


# def PlotMoments(fnl0Object:Histogram, fnl100Object:Histogram):
#     # fA = open("./hist/{f}_moments.txt".format(f=fnlAsimx), "r")
#     # fB = open("./hist/{f}_moments.txt".format(f=fnlBsimx), "r")

#     if not os.path.exists("./hist_plots"):    
#         os.makedirs("./hist_plots")
#     if not os.path.exists("./hist"):    
#         os.makedirs("./hist")

#     fnl0Object.calc_moment()
#     fnl100Object.calc_moment()

#     s = np.linspace(fnl0Object.sMin, fnl0Object.sMax, fnl0Object.sBins) #  TO DO: is it sBins or sBins-1?

#     # MOMENT 1
#     fig1, axs1 = plt.subplots(2)
#     axs1[0].plot(s, fnl0Object.m1, label="fnl 0")
#     axs1[0].plot(s, fnl100Object.m1, label="fnl 100")

#     axs1[1].plot(s, fnl100Object.m1 - fnl0Object.m1)

#     fig1.suptitle("{d} and {r} moment 1".format(d=fnl0Object.data_file, r=fnl100Object.data_file))
#     fig1.xlabel("s" +  r" $h^{-1}$Mpc")
#     axs1[0].ylabel("W0W1/<B0B1> moments")
#     axs1[1].ylabel("fnl100 moment - fnl0 moment")          
#     fig1.legend()
#     fig1.savefig("./hist_plots/sim{f}_moment1.png".format(f=fnl0Object.data_file[7:]))

#     # MOMENT 2
#     fig2, axs2 = plt.subplots(2)
#     axs2[0].plot(s, fnl0Object.m2, label="fnl 0")
#     axs2[0].plot(s, fnl100Object.m2, label="fnl 100")

#     axs2[1].plot(s, fnl100Object.m2 - fnl0Object.m2)

#     fig2.suptitle("{d} and {r} moment 1".format(d=fnl0Object.data_file, r=fnl100Object.data_file))
#     fig2.xlabel("s" +  r" $h^{-1}$Mpc")
#     axs2[0].ylabel("W0W1/<B0B1> moments")
#     axs2[1].ylabel("fnl100 moment - fnl0 moment")          
#     fig2.legend()
#     fig2.savefig("./hist_plots/sim{f}_moment2.png".format(f=fnl0Object.data_file[7:]))

#     # MOMENT 3
#     fig3, axs3 = plt.subplots(2)
#     axs3[0].plot(s, fnl0Object.m3, label="fnl 0")
#     axs3[0].plot(s, fnl100Object.m3, label="fnl 100")

#     axs3[1].plot(s, fnl100Object.m3 - fnl0Object.m3)

#     fig3.suptitle("{d} and {r} moment 1".format(d=fnl0Object.data_file, r=fnl100Object.data_file))
#     fig3.xlabel("s" +  r" $h^{-1}$Mpc")
#     axs3[0].ylabel("W0W1/<B0B1> moments")
#     axs3[1].ylabel("fnl100 moment - fnl0 moment")        
#     fig3.legend()
#     fig3.savefig("./hist_plots/sim{f}_moment3.png".format(f=fnl0Object.data_file[7:]))

#     # MOMENT 4
#     fig4, axs4 = plt.subplots(2)
#     axs4[0].plot(s, fnl0Object.m4, label="fnl 0")
#     axs4[0].plot(s, fnl100Object.m4, label="fnl 100")

#     axs4[1].plot(s, fnl100Object.m4 - fnl0Object.m4)

#     fig4.suptitle("{d} and {r} moment 1".format(d=fnl0Object.data_file, r=fnl100Object.data_file))
#     fig4.xlabel("s" +  r" $h^{-1}$Mpc")
#     axs4[0].ylabel("W0W1/<B0B1> moments")
#     axs4[1].ylabel("fnl100 moment - fnl0 moment")        
#     fig4.legend()
#     fig4.savefig("./hist_plots/sim{f}_moment4.png".format(f=fnl0Object.data_file[7:]))


    # plt.figure(2)
    # plt.plot(s, fnl0Object.m2, label="fnl 0")
    # plt.plot(s, fnl100Object.m2, label="fnl 100")

    # plt.title("{d} and {r} moment 2".format(d=fnl0Object.data_file, r=fnl100Object.data_file))
    # plt.xlabel("s" +  r" $h^{-1}$Mpc")
    # plt.ylabel("W0W1/<B0B1> moments log scale")        
    # plt.yscale("log")
    # plt.legend()
    # plt.savefig("./hist_plots/sim{f}_moment2.png".format(f=fnl0Object.data_file[7:]))

    # plt.figure(3)
    # plt.plot(s, fnl0Object.m3, label="fnl 0")
    # plt.plot(s, fnl100Object.m3, label="fnl 100")

    # plt.title("{d} and {r} moment 3".format(d=fnl0Object.data_file, r=fnl100Object.data_file))
    # plt.xlabel("s" +  r" $h^{-1}$Mpc")
    # plt.ylabel("W0W1/<B0B1> moments log scale")        
    # plt.yscale("log")
    # plt.legend()
    # plt.savefig("./hist_plots/sim{f}_moment3.png".format(f=fnl0Object.data_file[7:]))

    # plt.figure(4)
    # plt.plot(s, fnl0Object.m4, label="fnl 0")
    # plt.plot(s, fnl100Object.m4, label="fnl 100")

    # plt.title("{d} and {r} moment 4".format(d=fnl0Object.data_file, r=fnl100Object.data_file))
    # plt.xlabel("s" +  r" $h^{-1}$Mpc")
    # plt.ylabel("W0W1/<B0B1> moments log scale")        
    # plt.yscale("log")
    # plt.legend()
    # plt.savefig("./hist_plots/sim{f}_moment4.png".format(f=fnl0Object.data_file[7:]))


# if __name__ == "__main__":
#     ProcessHistogram('fnl0sim1.fits', 'randoms.fits', 'ex_par.txt', 45)
