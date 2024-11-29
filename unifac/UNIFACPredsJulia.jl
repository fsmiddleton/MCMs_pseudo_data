# use the below commands to install the necessary packages 
import Pkg; Pkg.add("TensorCast")
using Pkg; Pkg.add("XLSX")
import Pkg; Pkg.add("Clapeyron")
# Import packages 
using Clapeyron
using TensorCast
import XLSX

#Import compound names to be used 
xf = XLSX.readxlsx("unifac/UNIFACParams.xlsx")
compoundnames = xf["compounds"]["C1:C108"] 


#Create predictions and export them for all temperature and compound combinations 
lengthX = 19
x = range(0,1,length=lengthX)
X = Clapeyron.FractionVector.(x)
Temps = [293.15; 298.15; 303.15; 308.15; 313.15; 318.15; 323.15]
for T in Temps
    println(T)
    mixture = ["compound1" "compound2"]
    HE = []
    for comp1 in compoundnames 
        for comp2 in compoundnames
            #find whether it is a true mixture 
            if (cmp(comp1,comp2)!=0) 
                try
                    model = UNIFAC([comp1, comp2]) #create model
                    mixture = [mixture; comp1 comp2] #save mixture
                    append!(HE, [mixing.(model,1.013e6,T,X,enthalpy)]) #save prediction of excess enthalpy
                    println("Mixture predicted using UNIFAC: ",comp1, "+", comp2)
                catch 
                    println("Mixture did not exist in parameters for UNIFAC: ",comp1, "+", comp2)
                    try 
                        model = PR([comp1, comp2]) #create model
                        mixture = [mixture; comp1 comp2] #save mixture
                        append!(HE, [mixing.(model,1.013e6,T,X,enthalpy)])
                        println("Mixture predicted using PR: ",comp1, "+", comp2)
                    catch 
                        println("Mixture not in UNIFAC or PR parameters: ",comp1, "+", comp2)
                    end 
                end
            end 
        end 
    end 
    #export data 
    #create data as table
    columns = Vector()
    push!(columns,mixture[2:end,1])
    push!(columns,mixture[2:end,2])
    #first line of the 2d array 
    temp = [HE[1,1][1]]
    for m = 2:lengthX
        temp = [temp HE[1,1][m]]
    end 
    Hetemp = temp
    #rest of the array of HE data
    for j = 2:size(HE,1) #number of mixtures
        temp = [HE[j,1][1]]
        for m = 2:1lengthX
            temp = [temp HE[j,1][m]]
        end 
        Hetemp = [Hetemp; temp]
    end 
    #finish creating columns
    for k = 1:size(Hetemp,2)
        push!(columns, Hetemp[:,k])
    end
    #labels for excel sheet 
    labels = Any[1 2 x']
    #write data to file 
    XLSX.openxlsx(string("data/UNIFACPreds",string(T),".xlsx"), mode="w") do xf
        sheet = xf[1]
        XLSX.rename!(sheet, string(T))
        XLSX.writetable!(sheet,columns,labels)
    end

    println("Exported")
end 