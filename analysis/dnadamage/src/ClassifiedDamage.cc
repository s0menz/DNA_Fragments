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
	
	//declaration of vector for DNAFragments and map for Fragments size distribution
	std::vector<int> DNAFragments;	
	int binWidth = 200000;     //bin width set at 200 kBp
	std::map<int, int> FragmentSizeFreq;
	
	// ------------
	//declaration of vector for Cluster Sizes  and map for dsbClusters size distribution
	std::vector<int> ClusterSizes;
	int binWidth2 = 1;        //bin width set at 1 bp
	std::map<int, int> ClusterSizeFreq;
	
	
	// ----------
	
	for(int i = 0; i < dsbClusters.size(); i++){
		auto& a = dsbClusters[i];
		auto& b = dsbClusters[i+1];
		
		if(a.GetClusterChromo() == b.GetClusterChromo()){
		
		//DNA fragments	size
		int distance = b.GetBpBegin() - a.GetBpEnd();
		//Dsb Clusters size 	
		int ClusterSize = a.GetBpEnd() - a.GetBpBegin();
		
		DNAFragments.push_back(distance);
		ClusterSizes.push_back(ClusterSize);
	    
	    
	    //write the sorted list of fragments
	   	std::ofstream file("/home/local1/build-dsbandrepair_dev/_cluster_distances.txt", std::ios::out | std::ios::app);  
		file << "Distance between cluster " << i
            << " (begin=" << a.GetBpBegin() << ", " 
            << " end=" << a.GetBpEnd() << ") and cluster " << i+1
            << " (begin=" << b.GetBpBegin() << ", "
            << " end="    << b.GetBpEnd() << "): "
            << distance << " bp" << ", " 
            << "in chromosome: " << a.GetClusterChromo() << std::endl;

		   }
        
            
        }  
	
	
	
	// fill the map for fragments sizes distribution
     for (auto& distance : DNAFragments) {
		
		int bin = (distance / binWidth) * binWidth;
		FragmentSizeFreq[bin]++;	
			
	           }
	           
	 // fill the map for dsb Cluster sizes distribution
	  for (auto& ClusterSize : ClusterSizes) {
		  
		  int bin  = (ClusterSize / binWidth2) * binWidth2;
		  ClusterSizeFreq[bin]++;
	  } 
	  
	         
     //write the distribution of Fragment sizes 
        std::ofstream distrib_file("/home/local1/build-dsbandrepair_dev/_Fragments_distribution.txt", std::ios::out | std::ios::app);      
        for (auto& pair : FragmentSizeFreq){
			
               distrib_file << pair.first << " - " 
                            << pair.first + binWidth -1 << "	"
                            << pair.second << std::endl;  
       }
       
       //write the distribution of Cluster sizes 
        std::ofstream distrib_file2("/home/local1/build-dsbandrepair_dev/_Cluster_distribution.txt", std::ios::out | std::ios::app);      
        for (auto& pair : ClusterSizeFreq){
			
               distrib_file2 << pair.first << " - " 
                            << pair.first + binWidth2 -1 << "	"
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
