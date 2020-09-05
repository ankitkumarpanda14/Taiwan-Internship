# Import the libraries for use
import pandas as pd
import numpy as np
from datetime import datetime
from statistics import mean

"""
Variable Names and Definition:
    data: For reading the complete data from the raw data csv file
    hv: Series variable for storing the hot valve usage data
    wv: Series variable for storing the warm valve usage data
    cv: Series variable for storing the cold valve usage data
    VolumeUse: Dataframe for saving the volume used and the output csv is obtained from this dataframe
    dheatloss: Dataframe for getting the heat loss coefficient 
    heatlossindex: list variable to keep the index values of the dataframe dheatloss
    heatlossindcont: list variable to keep the continuous index values of the dataframe dheatloss
    heati: Selecting the continuous array i from heatlossindcont
    m: Series variable for storing the slope of each array
    n: The average of the series m giving the conductivity value
    dcoldloss: Dataframe for getting the cold loss coefficient 
    coldlossindex: list variable to keep the index values of the dataframe dcoldloss
    coldlossindcont: list variable to keep the continuous index values of the dataframe coldloss
    coldi: Selecting the continuous array i from coldlossindcont
    m1: Series variable for storing the slope of each array
    n1: The average of the series m giving the conductivity value
    dfwarmloss: Dataframe for getting the warm tank conductivity 
    warmlossindex: list variable to keep the index values of the dataframe dfwarmloss
    warmlossindcont: list variable to keep the continuous index values of the dataframe dfwarmloss
    warmi: Selecting the continuous array i from heatlossindcont
    m1: Series variable for storing the slope of each array
    n2: The average of the series m giving the conductivity value
    dheat: Data frame for determining the minimum and maximum heat temperatures for regular heating
    ascendingorderObjhot: dheat temperature values sorted in ascending order
    htcount: Counts of each temperature of the dheat dataframe
    heatmax: Upper limit temperature
    heatmin: Lower limit temperature
    dcool: Data frame for determining the minimum and maximum heat temperatures for regular heating
    ascendingorderObjc: dheat temperature values sorted in ascending order
    ctcount: Counts of each temperature of the dcool dataframe
    coldmax: Upper limit temperature
    coldmin: Lower limit temperature
"""

# Import and read data, drop the fields not needed for parameter calculation
""" Change the file name here for other simulators """

data = pd.read_csv("MA_04_01_test (1).csv", header = 0) 

data = data.drop([ 'Watts',  'HighWaterLevel', 'LowWaterLevel', 
                  'MeanWaterLevel', 'Sterilizing', 'ErrorCode'], axis = 1)



## Change the upload times to Unix Timestamps
def DT(row):
    date = datetime.strptime(row, "%Y-%m-%d %H:%M")
    return datetime.timestamp(date)
data['UploadTime'] = data['UploadTime'].apply(DT)

hv = (data['Hot_Valve'] > 0)  # To get the rows where the hot tank water is used
wv = (data['Warm_Valve'] > 0) # To get the rows where the warm tank water is used
cv = (data['Cold_Valve'] > 0) # To get the rows where the cold tank water is used

# Volume used Calculation
VolumeUse = pd.DataFrame(columns = ['VolumeUse', 'Case'])

for i, row in data.iterrows(): # Calculating the volume used by subtracting subsequent rows
    VolumeUse['VolumeUse'] = (-(data['Usage_L'] + 0.001*data['Usage_CC']) + (data['Usage_L'].shift(-1) + 0.001*(data['Usage_CC'].shift(-1))))

for j in range(len(hv)): # Setting the number as per the tank from which the water is used
    if(hv[j] == True):
        VolumeUse['Case'][j] = 1
    elif(wv[j] == True):
        VolumeUse['Case'][j] = 2
    elif(cv[j] == True):
        VolumeUse['Case'][j] = 3
    else:
        VolumeUse['Case'][j] = 0        

VolumeUse.Case = VolumeUse.Case.shift(-1)

for i,row in VolumeUse.iterrows(): # Setting the rows where volume use exist but the tank information is not present 
    if((VolumeUse['VolumeUse'][i] > 0) & (VolumeUse['Case'][i] == 0)):
        VolumeUse['Case'][i] = 1

# Generating a .csv file for the volume of water used and the case for the cpp program
""" Change the name of the file and path here """

VolumeUse.to_csv(r'C:\Users\AnkitPanda\Desktop\Internship@Taiwan\Test1.csv', index = None, header = None) 

# Heat Loss Coefficient for the hot tank 

# Getting the rows having the Heating flag as 0
dheatloss = pd.DataFrame()
for i, row in data.iterrows(): 
    if(row['Heating'] == 0):
        dheatloss = dheatloss.append(data.loc[i])

# Dropping coloumns not needed for calculation
dheatloss = dheatloss.drop(['ColdTemp', 'Cooling', 'Refilling', 'WarmTemp'], axis = 1)

# Separating the rows where the heat loss occurs 
heatlossindex = dheatloss.index[dheatloss['Heating'] == 0].tolist() # The index of the main dataframe where the heating is 0
heatlossindcont = []
k = 0
li = []
length = len(heatlossindex)
for i,item in enumerate(heatlossindex): # Separating out the continuous indexes to a new list element 
    if(i+1 < length and item + 1 != heatlossindex[i+1]):
        li.append(item)
        heatlossindcont.append(li)
        li = []
        continue
    li.append(item)

# Removing the rows which have the heating after refilling
heatarray = []
for i in range(len(heatlossindcont)):
    if(len(heatlossindcont[i]) >= 30): # Removing the list elements in which the refilling occurs, 30 mins. has been set as a time for the threshold 
        heatarray.append(heatlossindcont[i]) # If no one use the dispenser for more than 30 mins. only that data to be used

heatlossindcont = heatarray
        
# Fitting a linear line to the data for obtaining the conductivity value
m = []
for i in range(len(heatlossindcont)):
    heati = heatlossindcont[i]
    l = len(heatlossindcont[i])
    f = heati[0]
    h = heati[l - 1]
    a = data.iloc[f:h,1]
    b = list(range(len(a)))
    c = np.polynomial.polynomial.Polynomial.fit(b,a,1)
    m.append(c.convert().coef[1]) # Getting the slope of the lines after fitting the curve

n = mean(m) # Calculating the mean of the slope data to get the conductivity


# Getting the rows where the value of cooling is 0 
dcoldloss = pd.DataFrame()
for i, row in data.iterrows(): 
    if(row['Cooling'] == 0):
            dcoldloss = dcoldloss.append(data.loc[i])

# Dropping the rows which do not have use in the calculation 
dcoldloss = dcoldloss.drop(['Heating', 'HotTemp', 'Refilling', 'WarmTemp'], axis = 1)

# Getting the list of index which have cooling flag is 0
coldlossindex = dcoldloss.index[dcoldloss['Cooling'] == 0].tolist()
coldlossindcont = []
li1 = []
length = len(coldlossindex)
for i,item in enumerate(coldlossindex): # Separating the continuous rows
    if(i+1 < length and item + 1 != coldlossindex[i+1]):
        li1.append(item)
        coldlossindcont.append(li1)
        li1 = []
        continue
    li1.append(item)

# Removing the coloumns which have cold water usage
coldarray = []
for i in range(len(coldlossindcont)):
    if(len(coldlossindcont[i])> 100): # If a length of less than 100 minutes means tha the cold water was used and hence cannot be tested 
        coldarray.append(coldlossindcont[i])
        
coldlossindcont = coldarray
    
m1 = []

# Fitting a linear line to the points for the conductivity value
coldi = coldlossindcont[0]
l1 = len(coldlossindcont[0])
f1 = coldi[0]    
h1 = coldi[l1 - 1]
a1 = data.iloc[f1:h1+5,10]
b1 = list(range(len(a1)))
c1 = np.polynomial.polynomial.Polynomial.fit(b1,a1,1)
m1.append(c1.convert().coef[1])

n1 = mean(m1)

" Warm tank coefficient "
 
dfwarmloss = pd.DataFrame() # Dataframe for warm temperature values 
for i, row in data.iterrows(): # getting the values where the refilling is not happening
    if((row['Refilling'] == 0) & (row['Warm_Valve'] == 0)):
        dfwarmloss = dfwarmloss.append(data.loc[i])


warmlossindex = dfwarmloss.index[(dfwarmloss['Refilling'] == 0) & (dfwarmloss['Warm_Valve'] == 0)].tolist() # getting the index values of the warm tank

dfwarmloss = dfwarmloss.drop(['ColdTemp', 'Cooling', 'Heating', 'HotTemp'], axis = 1)

warmlossindcont = []
li1 = []
length = len(warmlossindex)
for i,item in enumerate(warmlossindex): # Getting the row indices into a new list
    if(i+1 < length and item + 1 != warmlossindex[i+1]):
        li1.append(item)
        warmlossindcont.append(li1)
        li1 = []
        continue
    li1.append(item)

warmarray = []
for i in range(len(warmlossindcont)): # If the warm tank is refilling in less than 20 minutes then we need to remove those rows
    if(len(warmlossindcont[i]) > 20):
        warmarray.append(warmlossindcont[i])
        
warmlossindcont = warmarray

m2 = []
for i in range(len(warmlossindcont)):
    warmi = warmlossindcont[i]
    l1 = len(warmlossindcont[i])
    f1 = warmi[0]
    h1 = warmi[l1 - 1]
    a1 = data.iloc[f1:h1,9]
    b1 = list(range(len(a1)))
    c1 = np.polynomial.polynomial.Polynomial.fit(b1,a1,1)
    m2.append(c1.convert().coef[1])

mavg = []
for i in range(len(m2)): # Removing the lines which have a negative slope for the equation
    if(m2[i] < 0):
        mavg.append(m2[i])
        
n2  = mean(mavg) # Calculating the mean of all the conductivity values

    
" Getting the insulation temperature for the hot tank "

# Separating the columns in which the Heating, Power Saving and Refilling is not there
dheat = pd.DataFrame()
for i, row in data.iterrows():
    if((row['Heating'] == 0) & (row['SavingPower'] == 0) & (row['Refilling'] == 0)):
        dheat = dheat.append(data.loc[i])

# Sorting the hot temperatures to get the minima and maxima of the data
ascendingorderObjhot = dheat.sort_values(by = 'HotTemp')

# Dropping the unwanted rows
ascendingorderObjhot = ascendingorderObjhot.drop(['ColdTemp', 'Cooling', 'Heating', 'Hot_Valve', 
                                            'Refilling', 'SavingPower', 'UploadTime', 'WarmTemp', 'Warm_Valve' ], axis = 1)

# Getting the value counts of the objects     
htcount = ascendingorderObjhot['HotTemp'].value_counts()

# Sorting the index by temperatures
htcount = htcount.sort_index()

# Removing the values which are not frequent in the data
# The frequency here has been set to 10 since observed in the dataset the extreme values have a frequency of 4 or 5
htcount = htcount[htcount >= 10]

# Getting the minima and maxima of the data
heat = []
heat = htcount.index
heatmin = heat[0]
l = len(heat)
heatmax = heat[l - 1]

" Getting the insulation temperature for the cold tank"

# Separating the columns in which the Cooling, PowerSaving and Cold_Valve is not used
dcool = pd.DataFrame()
for i,row in data.iterrows():
    if((row['Cooling'] == 0) & (row['SavingPower'] == 0) & (row['Cold_Valve'] == 0)):
        dcool = dcool.append(data.loc[i])

# Sorting the hot temperatures to get the minima and maxima of the data
ascendingorderObjc = dcool.sort_values(by = 'ColdTemp')

# Dropping the unwanted columns
ascendingorderObjc = ascendingorderObjc.drop(['HotTemp', 'Cooling', 'Heating', 'Hot_Valve', 
                                            'Refilling', 'SavingPower', 'UploadTime', 'WarmTemp', 'Warm_Valve', 'Cold_Valve' ], axis = 1)

# Getting the value counts of the objects     
ctcount = ascendingorderObjc['ColdTemp'].value_counts()

# Sorting the index by temperatures
ctcount = ctcount.sort_index()

# Removing the values which are not frequent in the data
# The frequency here has been set to 10 since observed in the dataset the extreme values have a frequency of 4 or 5
ctcount = ctcount[ctcount >= 10]

# Getting the minima and maxima of the data
cold = []
cold = ctcount.index
coldmin = cold[0]
lc = len(cold)
coldmax = cold[lc - 1]


" Write the final parameters to the csv for the simulator program to read it "
par = pd.DataFrame({'0' : [heatmin], '1' : [heatmax], '2': [coldmin], '3': [coldmax], '4' : [n],'5' : [n2], '6' : [n1]})
par.to_csv(r'C:\Users\AnkitPanda\Desktop\Internship@Taiwan\parameter.csv', index = None, header = None)
