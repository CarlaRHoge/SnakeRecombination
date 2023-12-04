%matplotlib inline
f, (a0) = plt.subplots(1, figsize=(7,6))
points = pd.DataFrame()
#For each "Mean" output of Find_Means_Boots_Counts.py for this plot (make sure to update the color "c" and label "l")
for file, c, l, in zip([os.path.join(folder, "LDHelmet", "Features", "PRDM9_All_Top50kUpd_all_Means_dict.txt"), os.path.join(folder, "LDHelmet", "Features", "Promoters_all_Means_dict.txt")], ["blue", "orange"], ["PRDM9 binding sites", "Promoter-like features"]): 
    feature = file.split("/")[-1].split("_")[0]
    boot = file.split("Mean")[0] + "Boot" + file.split("Mean")[-1]                
    with open(boot, "rb") as f:
        boots = pickle.load(f)
    #for every bootstrap iteration, find a lowess curve. Note that the 101 is the 1+number of x values in the final plot - in this case, 100 bins up to 10000 is 100 bins. 
    smoothed_values = np.empty((len(boots[0]), 101))
    for count in range(len(boots[0])): 
        x = []
        y = []
        for dis in sorted(boots): 
            if dis > 10000: continue
            x.append(dis)
            y.append(boots[dis][count])
        #Smooth the bootstrapped values for this run
        smoothed_values[count] = lowess(y, x, frac=.1, it=3, return_sorted = False)
    #Sort each column
    sorted_values = np.sort(smoothed_values, axis=0)
    #Find the confidence intervals for each distance
    bound = int(len(boots[0]) * (1 - .95) / 2)
    bottom = sorted_values[bound-1]
    top = sorted_values[-bound]
    x=[]
    y=[]
    y_c=[]
    x_c=[]
    lows = 0
    clows = 0
    highs = 0
    chighs = 0
    with open(file, "rb") as f:
        means = pickle.load(f)
    #For each distance in the "means" file
    for dis in sorted(means): 
        #Skip it if it's too large
        if dis > 10000: continue
        #Work towards finding the average if it's >8000bp, to make the plot relative to rec rate at 8-10kb
        if dis > 8000: 
            highs += means[dis]
            chighs += 1
        #add distance to x and mean rec rate to y.
        x.append(dis)
        y.append(means[dis])
    #The following code makes a histogram of the number of values at each distance, if you want it. 
    cou = file.split("Mean")[0] + "Count" + file.split("Mean")[-1]     
    with open(cou, "rb") as f:
        counts = pickle.load(f)
    newcounts = {}
    for dis in counts: 
        if dis > 10000: continue
        newdis = dis//20 * 20
        if newdis in newcounts: 
            newcounts[newdis] += counts[dis]
        else: 
            newcounts[newdis] = counts[dis]
    for dis in newcounts: 
        x_c.append(dis)
        y_c.append(newcounts[dis])
    #Find average rec rate at 8-10kb 
    highs = highs / chighs
    #Divide each y value by average rec rate at 8-10kb 
    y = [i/highs for i in y]
    bottom = [i/highs for i in bottom]
    top = [i/highs for i in top]
    #Smoooth
    lw = lowess(y, x, frac=.1, it=3, return_sorted = False)
    #Plot confidence interval
    a0.fill_between(x, bottom, top, alpha = .1, color=c)
    #Plot actual curve 
    a0.plot(x, lw, linewidth=2, color=c, label=l)
    a0.set_ylabel("Relative Recombination Rate", fontsize=14)
    a0.set_xlabel("Distance to Feature (bp)", fontsize=14)
    a0.legend(fontsize=13)
#a0.set_ylim((.75, 4.5))
plt.title("All alleles No CpGI/TSS", fontsize=20)
plt.savefig(os.path.join(folder, "LDHelmet", "Features", "RecVsDistance_PRDM9All-Promoters_10kb.png"))