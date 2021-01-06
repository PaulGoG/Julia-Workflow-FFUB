using DataFrames, DataFramesMeta
using CSV

#reading the CSV file
cd("C:\\Users\\GoG\\Desktop\\JuliaCodes\\Dubna2019Cosmology\\Craciun")
df = CSV.read("data\\RawData.csv")

#taking the columns that are needed
ex = df[12:end,[8,5,11]]
#dropping df to clear memory
df = nothing
#formatting the dataframes obj as needed
colnames = [ex[1,1],ex[1,2],"redshift"]
names!(ex,Symbol.(colnames))
ex = ex[2:end,:]
#removing rows with missing data
dropmissing!(ex)

#selecting type Ia supernovae
ex = @linq ex|>where("SNIa".==:Method)
#deleting duplicate data
ex = @linq ex|>where(nonunique(ex,3).==false)
#sorting by redshift
sort!(ex,cols=[:redshift])
#removing the method column
ex = ex[:,[2,3]]

#writing a new CSV file
CSV.write("data\\prelucratedData.csv",ex)
