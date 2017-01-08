# DNA per nucleus 
dsDNAwtPerNuc <- 607.4 #g/mol or Da
Avogadro <- 6.022e+23

NtInCeGenome <-1e+8 #1 Mb DNA in genome
 
g2ug<-1e+6 # number of ug in a gram

# Mwt of Ce genome in g/mole
ceGenome_gPerMole <- NtInCeGenome*dsDNAwtPerNuc
ceGenome_gPerMole

# Mwt of a single genome (Mwt/ avogaro's number)
singleGenome_g <- ceGenome_gPerMole/Avogadro
singleGenome_g

#convert from g to ug
singleGenome_ug<- singleGenome_g*g2ug
singleGenome_ug

# genome is diploid in each nucleus -> mutliply by 2
diploidGenome_ug<-2*singleGenome_ug
diploidGenome_ug

# how many diploid genomes (nuclei) do you need for 5 ug of DNA?
NumNuc<-5/diploidGenome_ug
NumNuc

# express in millions of nuclei
NumMillionNuc<-NumNuc/1e+6
NumMillionNuc



########################
# 1 ug of 1000 bp DNA = 9.1x 10^11 molecules
# 1 ug of 1 Mb DNA = 9.1 x 10^8 molecules
# 1 ug of 100 Mb DNA (ceGenome) = 9.1x 10^6 molecules
# 5 ug of 100 Mb DNA = 4.55 x 10^7 molecules  i.e. 45 million

Wt1Mb<- 1/(9.1*10^8)
CeGenomeWt<- 100 *Wt1Mb*2
CeGenomeWt
5/CeGenomeWt/1e6
#23 million nuclei for 5 ug

HumanGenomeWt<-3000*Wt1Mb*2
HumanGenomeWt
6/HumanGenomeWt/1e6
#910,000 nuclei for 5 ug

flyGenomeWt<-180*Wt1Mb*2
flyGenomeWt
6/flyGenomeWt/1e6
# 15 million nuclei for 5 ug

HumanGenomeWt*250000

flyGenomeWt*2.5e6

# in 250,000 cells
cellNum=250000
cellNum*HumanGenomeWt

