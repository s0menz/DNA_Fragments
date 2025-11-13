//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
//
/// \file ClassifiedDamage.cc
/// \brief Implementation of the ClassifiedDamage class

#include "ClassifiedDamage.hh"
#include <fstream>
#include <iostream>
#include <vector>
#include <algorithm> 
#include <map>
//#include "G4SystemOfUnits.hh"


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

ClassifiedDamage::ClassifiedDamage()
{
	fType = fNA;
	fDamage.clear();
	fBp_begin = 0;
	fBp_end = 0;
	fComplexity = -1;
	fIncludeBase = false;
	
	
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void ClassifiedDamage::AddDamage(Damage pDmg)
{
	fDamage.push_back(pDmg);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void ClassifiedDamage::ComputeBp()
{
	fBp_begin = fDamage[0].GetCopyNb();
	fBp_end = fDamage[fDamage.size()-1].GetCopyNb();

	fBp_barycenter=0;
	for(auto it=fDamage.begin();it!=fDamage.end();it++)
	{
		fBp_barycenter+=it->GetCopyNb();
	}
	fBp_barycenter = fBp_barycenter/fDamage.size();
	fChromo = fDamage[0].GetChromo();
	
	
	  
 
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void ClassifiedDamage::ComputePosition()
{

}

std::vector<ClassifiedDamage> ClassifiedDamage::FragmentSize(const std::vector<ClassifiedDamage>& listClassifiedDamage) 
{  
	
   
   //new vector with only the dsb clusters
   std::vector<ClassifiedDamage> dsbClusters;
  


   
   for (auto& ClassDamage : listClassifiedDamage){
	     
	     //I work on a copy of the classified damages
         ClassifiedDamage copy = ClassDamage; 
	    if (copy.GetClassifiedDamageType() == ClassifiedDamage::ClassifiedDamageType::fDSB){
			copy.ComputeBp();
			dsbClusters.push_back(copy);
			
		}
		
  }		
  
  
  
  
				
 
return dsbClusters;
}		
   
     
     
     	
 



void ClassifiedDamage::SortDsbClusters(std::vector<ClassifiedDamage>& dsbClusters)
{
	//chromosomes lenghts for Fibroblast nucleus. 
	//Must be changed if you use an other DNAFabric nucleus
	std::vector<long double> chromoLenght = {  
			2.586260E+08, 2.585990E+08,	
		2.482510E+08,	2.483350E+08,	
		2.019270E+08,	2.020040E+08,	
		1.911760E+08,	1.911110E+08,	
		1.852220E+08,	1.849900E+08,	
		1.745840E+08,	1.744460E+08,	
		1.621350E+08,	1.621600E+08,	
		1.477830E+08,	1.478370E+08,	
		1.425710E+08,	1.427640E+08,	
		1.390800E+08,	1.392300E+08,	
		1.385650E+08,	1.384350E+08,	
		1.372760E+08,	1.370910E+08,	
		1.163270E+08,	1.178180E+08,	
		1.113470E+08,	1.113770E+08,	
		1.080560E+08,	1.083000E+08,	
		9.505560E+07,	9.504480E+07,	
		8.721980E+07,	8.743680E+07,	
		8.162030E+07,	8.167270E+07,	
		6.107600E+07,	6.102330E+07,	
		6.746970E+07,	6.760940E+07,	
		4.893800E+07,	4.900600E+07,	
		5.352630E+07,	5.370580E+07,	
		1.579820E+08,	1.578690E+08}; //the 45th and 46th chromosome are XX in this case
// ----------------------------------------------------------------------
// ----------------------------------------------------------------------
// --------------------  conversion base pairs --> nanometers -----------
// ----------------------------------------------------------------------
	
	double bpDist = 1;             //bp ---> nanometers =  0.34 nm/bp
// ----------------------------------------------------------------------
// ----------------------------------------------------------------------
// ----------------------------------------------------------------------
// ----------------------------------------------------------------------
	
	
    for (auto& ClassDamage : dsbClusters){
	      ClassDamage.ComputeBp();
	  }
	      //sort for number of chromosomes and BpBegin 
	      //important:  the copyNb of Bp is reset starting a new chromosome
	
    
		  std::sort(dsbClusters.begin(), dsbClusters.end(), 
             [](const ClassifiedDamage& a, const ClassifiedDamage& b){
				 if(a.GetClusterChromo() == b.GetClusterChromo()){
				 return a.GetBpBegin() < b.GetBpBegin();
				 }
				 return a.GetClusterChromo() < b.GetClusterChromo();
				});
				
	// Put the clusters in different vectors according to the chromosome			
	std::vector<std::vector<ClassifiedDamage>> dsbClustInChromosomes(46);
	
	for (auto& ClassDamage : dsbClusters){
		int chromo = ClassDamage.GetClusterChromo();
		
		dsbClustInChromosomes[chromo].push_back(ClassDamage);
	    }
	  
	    
    // just a check    
	for(int i=0; i < 46 ; i++){
		for(auto& ClassDamage : dsbClustInChromosomes[i]){
    std::ofstream file("/home/local1/build-dsbandrepair_dev/_cluster_in_chromosomes.txt", std::ios::out | std::ios::app);  
		file << "chromo " << i << " check chromosome:" << ClassDamage.GetClusterChromo()
            << ", classDamage_begin " << ClassDamage.GetBpBegin()  
            << ",  end=" << ClassDamage.GetBpEnd()   << std::endl;
		   
        }
	}
        
         
      
            
            
	//declaration of vector for DNAFragments and map for Fragments size distribution
	std::vector<int> DNAFragments;	
	int binWidth = 200000;     //bin width for fragments set at 200 kBp
	std::map<int, int> FragmentSizeFreq;
	
	// ------------
	//declaration of vector for Cluster Sizes  and map for dsbClusters size distribution
	std::vector<int> ClusterSizes;
	int binWidth2 = 1;        //bin width for DSB clusters set at 1 bp
	std::map<int, int> ClusterSizeFreq;
	
	
	// ----------
	
	int distance = -1;
	int ClusterSize = -1;
	
	// Fill the DNAFragments vector with the fragments 
	// at the beginning of each chromosome
	for(int j = 0; j < 46; j++){
		
		auto& b = dsbClustInChromosomes[j][0];
		distance = b.GetBpBegin();
		std::ofstream file("/home/local1/build-dsbandrepair_dev/_cluster_distances.txt", std::ios::out | std::ios::app);  
		file << "chromo: " << j 
		<< " first damage at: " << b.GetBpBegin() 
		<< " distance: " << distance 
		 << std::endl;   
		 DNAFragments.push_back(distance); 
		 
	 }
	 
	 // Fill the DNAFragments vector with the fragments 
	// at the End of each chromosome
	 for(int j = 0; j < 46; j++){
		int lastDamage =  dsbClustInChromosomes[j].size() - 1;
		auto& a = dsbClustInChromosomes[j][lastDamage];
		distance = chromoLenght[j] - a.GetBpEnd();
		std::ofstream file("/home/local1/build-dsbandrepair_dev/_cluster_distances.txt", std::ios::out | std::ios::app);  
		file << "chromo: " << j 
		<< " last damage at: " << a.GetBpEnd() 
		<< " and chromo Lenght: " << chromoLenght[j] 
		<< " distance: " << distance 
		 << std::endl;   
		 DNAFragments.push_back(distance); 
	 }
	 
	 // Fill the DNAFragments vector and the ClusterSizes vector
	//  with fragm. and clusters in the middle of the chromosomes
	 for(int j = 0; j < 46; j++){
		for(int i = 0; i < dsbClustInChromosomes[j].size() - 1; i++){
		auto& a = dsbClustInChromosomes[j][i];
		auto& b = dsbClustInChromosomes[j][i+1];
		distance = b.GetBpBegin() - a.GetBpEnd();
		ClusterSize = a.GetBpEnd() - a.GetBpBegin();
		
		std::ofstream file("/home/local1/build-dsbandrepair_dev/_cluster_distances.txt", std::ios::out | std::ios::app);  
		file << "chromo: " << j 
		<< " damage " << i << " at: " << a.GetBpEnd() 
		<< " anddamage: " << i+1 << " at: " << b.GetBpBegin()
		<< " distance: " << distance 
		 << std::endl;   
		 DNAFragments.push_back(distance); 
		 ClusterSizes.push_back(ClusterSize);
	 }
 }
	/*for(int i = -1; i < dsbClustInChromosomes[j].size()+1; i++){
		
		
			if(i==-1){
				auto& b = dsbClustInChromosomes[j][0];
		
				distance = b.GetBpBegin();
			
				
	       } else if(i == dsbClustInChromosomes[j].size() ){
			   
                auto& a = dsbClustInChromosomes[j][i];
				distance = chromoLenght[j] - a.GetBpEnd();
                
               
		
		   } else if(i >-1 && i < dsbClustInChromosomes[j].size()  ){	 
			 auto& a = dsbClustInChromosomes[j][i]; 
			 auto& b = dsbClustInChromosomes[j][i+1];
		//DNA fragments	size
		distance = b.GetBpBegin() - a.GetBpEnd();
		
		//Dsb Clusters size 	
		ClusterSize = a.GetBpEnd() - a.GetBpBegin();
		
		
		ClusterSizes.push_back(ClusterSize);
	    }
	   
	    DNAFragments.push_back(distance); 
	    
	   /* if(a.GetClusterChromo() < b.GetClusterChromo()) {
			
			distance = b.GetBpBegin();
			ClusterSize = a.GetBpEnd() - a.GetBpBegin();
			DNAFragments.push_back(distance);
		ClusterSizes.push_back(ClusterSize);
			}*/
	    
	    /*
	    //Devi riscrivere questo check perche non ci sono piu a e b definiti!
	    //write the sorted list of fragments
	   	std::ofstream file("/home/local1/build-dsbandrepair_dev/_cluster_distances.txt", std::ios::out | std::ios::app);  
		file << "Distance between cluster " << i
            << " (begin=" << a.GetBpBegin() << ", " 
            << " end=" << a.GetBpEnd() << ") and cluster " << i+1
            << " (begin=" << b.GetBpBegin() << ", "
            << " end="    << b.GetBpEnd() << "): "
            << distance << " bp" << ", " 
            << "in chromosome: " << a.GetClusterChromo() << "first fragment: " << distanceFirst << " last fragment:" << distanceLast << std::endl;

        */
      /*    std::ofstream file("/home/local1/build-dsbandrepair_dev/_cluster_distances.txt", std::ios::out | std::ios::app);  
		file << "chromo: " << j 
		<< " i: " << i 
		<< " distance: " << distance 
		 << std::endl;   
        }  //closing the "i-for cicle"    
        
        */
	 // closing the "j-for cicle"
// -------------------------------------------------------------------
// -------------------------------------------------------------------
// ------------- maps for distributions     --------------------------
// -------------------------------------------------------------------
// -------------------------------------------------------------------
	
// Fill the map for Fragment sizes
    for (auto& distance : DNAFragments) {
		
		int bin = (distance / binWidth) * binWidth;
		FragmentSizeFreq[bin]++;		
	           }
	       
	           
 //adding the bins with 0 entries   
std::vector<int> binsToAddToFrag; // first, bins added to a vector

if (!FragmentSizeFreq.empty()) {
    auto it = FragmentSizeFreq.begin();
    auto prev = it;
    ++it;

    for (; it != FragmentSizeFreq.end(); ++it, ++prev) {
        int a = prev->first;
        int b = it->first;

        int gap = b - a;
        if (gap > binWidth) {
            for (int j = a + binWidth; j < b; j += binWidth) {
                binsToAddToFrag.push_back(j); // record the missing bin
            }
        }
    }

    // then, adding the bins to the map
    for (int bin : binsToAddToFrag) {
        FragmentSizeFreq[bin] = 0;
    }
}
 // --------------------------------------------------------------------
 // fill the map for dsb Cluster sizes distribution
	  for (auto& ClusterSize : ClusterSizes) {
		  
		  int bin  = (ClusterSize / binWidth2) * binWidth2;
		  ClusterSizeFreq[bin]++;
	  } 

        //adding the bins with 0 entries   
std::vector<int> binsToAddToCluster; // first, bins added to a vector

if (!ClusterSizeFreq.empty()) {
    auto it = ClusterSizeFreq.begin();
    auto prev = it;
    ++it;

    for (; it != ClusterSizeFreq.end(); ++it, ++prev) {
        int a = prev->first;
        int b = it->first;

        int gap = b - a;
        if (gap > binWidth2) {
            for (int j = a + binWidth2; j < b; j += binWidth2) {
                binsToAddToCluster.push_back(j); // memorizza il bin mancante
            }
        }
    }

    // then, adding the bins to the map
    for (int bin : binsToAddToCluster) {
        ClusterSizeFreq[bin] = 0;
    }
}
	
	 	   
     //write the distribution of Fragment sizes 
        std::ofstream distrib_file("/home/local1/build-dsbandrepair_dev/_Fragments_distribution.txt", std::ios::out | std::ios::app);      
        for (auto& pair : FragmentSizeFreq){
			
			   long double BinMinNM =  pair.first * bpDist;
			   long double BinMaxNM = (pair.first + binWidth -1) * bpDist;
			   
               distrib_file << BinMinNM << " - " 
                            << BinMaxNM << "		"
                            << pair.second << std::endl;  
       }
       
       //write the distribution of Cluster sizes 
        std::ofstream distrib_file2("/home/local1/build-dsbandrepair_dev/_Cluster_distribution.txt", std::ios::out | std::ios::app);      
        for (auto& pair : ClusterSizeFreq){
			   long double BinMinNM =  pair.first * bpDist;
			   long double BinMaxNM = (pair.first + binWidth2 -1) * bpDist;
			   
               distrib_file2 << BinMinNM << " - " 
                            << BinMaxNM << "	"
                            << pair.second << std::endl;  
       }
       
	std::ofstream logfile("/home/local1/build-dsbandrepair_dev/_dsbClusters.txt", std::ios::out | std::ios::app);	
	for (auto& a : dsbClusters){
		logfile  << "DSB cluster on chromosome: " 
		  << a.GetClusterChromo() << ","
          << a.GetBpBegin() << "," 
          << a.GetBpEnd() << "," 
          << " barycenter=" << a.GetBpBarycenter() 
          << std::endl;		}	 
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void ClassifiedDamage::ComputeComplexity()
{
	if(fDamage.size()==0)
	{
		fComplexity = -1;
	}
	else
	{
		fComplexity = -1;
		for(auto it=fDamage.begin();it!=fDamage.end();it++)
		{
			if((it->GetDamageType()==Damage::DamageType::fBackbone))
			{
				fComplexity++;
			}
			else
			{
				if(fIncludeBase)
				{
					fComplexity++;
				}
			}
		}
	}
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void ClassifiedDamage::ComputeType()
{
	bool firstStrandTouched = false;
	bool secondStrandTouched = false;

	if(fDamage.size()==0)
	{
		fType = fNA;
	}
	else
	{
		for(auto it=fDamage.begin();it!=fDamage.end();it++)
		{
			if((it->GetDamageType()==Damage::DamageType::fBackbone))
			{
				int strand = it->GetStrand();
				if(strand == 1)
				{
					firstStrandTouched = true;
				}
				if(strand == 2)
				{
					secondStrandTouched = true;
				}
			}

			if (it->GetCause() == Damage::DamageCause::fDirect) {
				fIsThereDirectContribution = true;
			}

			if (it->GetCause() == Damage::DamageCause::fIndirect) {
				fIsThereIndirectContribution = true;
			}
		}

		if(firstStrandTouched && secondStrandTouched)
		{
			fType = fDSB;
		}
		else
		{
			fType = fSSB;
		}
	}
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void ClassifiedDamage::Reset()
{
	fType = fNA;
	fDamage.clear();
	fBp_begin = 0;
	fBp_end = 0;
	fComplexity = -1;
}
