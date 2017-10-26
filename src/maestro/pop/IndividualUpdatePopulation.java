/**
 * Copyright 2017 Felipe Hernández
 * 
 * Licensed under the Apache License, Version 2.0 (the "License"); you may not use this file except
 * in compliance with the License. You may obtain a copy of the License at
 * 
 * http://www.apache.org/licenses/LICENSE-2.0
 * 
 * Unless required by applicable law or agreed to in writing, software distributed under the 
 * License is distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either
 * express or implied. See the License for the specific language governing permissions and
 * limitations under the License.
 */

package maestro.pop;

import java.util.ArrayList;
import java.util.HashSet;
import java.util.Iterator;

import probDist.Normal;
import utilities.stat.ContSeries;
import maestro.ContVar;
import maestro.DiscVar;
import maestro.solution.Solution;
import maestro.solution.SolutionWrapper;

/**
 * A population of candidate solutions that can be updated by offering a single solution at a time
 * @author Felipe Hernández
 */
public class IndividualUpdatePopulation implements Population 
{
	
	// --------------------------------------------------------------------------------------------
	// Constants
	// --------------------------------------------------------------------------------------------
	
	/**
	 * Default value for {@link #updateDist}
	 */
	public final static double D_UPDATE_DIST = 0.1;
	
	// --------------------------------------------------------------------------------------------
	// Attributes
	// --------------------------------------------------------------------------------------------
	
	/**
	 * Partitions in which solutions are stored. Each partition has a maximum size (number of 
	 * solutions) and a rarity standard (which all the solutions within should hold). Optionally, 
	 * the history of accepted solutions within the partition is kept.
	 */
	private ArrayList<Partition> partitions;
	
	/**
	 * The list of the discrete variables of the problem
	 */
	private ArrayList<DiscVar> discVars;
	
	/**
	 * The list of the continuous variables of the problem
	 */
	private ArrayList<ContVar> contVars;
	
	/**
	 * The solution capacity of the current population
	 */
	private int capacity;
	
	/**
	 * The number of new solutions refused from the population since the last solution was 
	 * accepted
	 */
	private int refuseCount;
	
	/**
	 * The number of processed solutions added to the population since the last update of the 
	 * distributions of the variables
	 */
	private int addedSinceUpdate;

	/**
	 * The percentage of new solutions accepted compared to the size of the population before the
	 * distributions of the variables are updated
	 */
	private double updateDist;
	
	// --------------------------------------------------------------------------------------------
	// Constructor
	// --------------------------------------------------------------------------------------------
	
	public IndividualUpdatePopulation(int capacity, int partitionCount, double minRarityStd, 
			double maxRarityStd, ArrayList<DiscVar> discVars, ArrayList<ContVar> contVars) 
	{
		this.capacity		= capacity;
		this.discVars		= discVars;
		this.contVars		= contVars;
		refuseCount			= 0;
		addedSinceUpdate	= 0;
		updateDist			= D_UPDATE_DIST;
		
		// Initialize partitions
		partitions			= new ArrayList<>();
		int accum 			= 0;
		int varCount 		= (discVars == null ? 0 : discVars.size()) 
								+ (contVars == null ? 0 : contVars.size());
		double adjustedMinRS = minRarityStd < maxRarityStd ? minRarityStd : maxRarityStd;
		double adjustedMaxRS = maxRarityStd > minRarityStd ? maxRarityStd : minRarityStd;
		if(varCount > 1)
		{
			Normal meanRarityDist = new Normal(0.5, 0.2897/Math.sqrt(varCount));
			adjustedMinRS = minRarityStd == 0.0 ? 0.0 : meanRarityDist.getInvCDF(minRarityStd);
			adjustedMaxRS = maxRarityStd == 1.0 ? 1.0 : meanRarityDist.getInvCDF(maxRarityStd);
			adjustedMinRS = adjustedMinRS < 0.0 ? 0.0 : adjustedMinRS;
			adjustedMaxRS = adjustedMaxRS > 1.0 ? 1.0 : adjustedMaxRS;
		}
		for(int i = 0 ; i < partitionCount ; i++)
		{
			int size = (int)(((double)(i + 1)/(double)partitionCount)*(double)capacity) - accum;
			size = i == 0 ? Math.max(2, size) : size;
			accum += size;
			double rarityStd = 0;
			if(i != 0)
			{
				if(partitionCount == 2)
					rarityStd = (adjustedMinRS + adjustedMaxRS)/2;
				else
					rarityStd = adjustedMinRS + ((double)(i - 1)/((double)(partitionCount - 2)))
													*(adjustedMaxRS - adjustedMinRS);
			}
			partitions.add(new Partition(size, rarityStd));
		}
	}
	
	// --------------------------------------------------------------------------------------------
	// Methods
	// --------------------------------------------------------------------------------------------

	@Override
	public int size() 
	{
		int size	= 0;
		for (Partition part : partitions)
			size	+= part.getSize();
		return size;
	}

	@Override
	public int getCapacity() 
	{
		return capacity;
	}
	
	@Override
	public void setCapacity(int capacity) 
	{
		this.capacity = capacity < 0 ? 0 : capacity;
	}
	
	@Override
	public double getUpdateDist() 
	{
		return updateDist;
	}

	@Override
	public void setUpdateDist(double updateDist) 
	{
		this.updateDist = updateDist < 0 ? 0 : updateDist;
	}
	
	@Override
	public Partition getPartition(int index) 
	{
		return partitions.get(index);
	}

	@Override
	public ArrayList<Partition> getPartitions() 
	{
		return partitions;
	}
	
	@Override
	public void forceAddSolution(SolutionWrapper solution) 
	{
		partitions.get(0).addSolution(solution);
		addedSinceUpdate++;
	}

	@Override
	public void offerSolution(SolutionWrapper solution) 
	{		
		// Check if the distributions of the variables should be updated
		if(solution.getGenIndex() != -1) // Verify if adding random solutions
			if(((double)addedSinceUpdate)/capacity >= Math.min(1.0, updateDist))
			{
				updateVarDist();
				addedSinceUpdate = 0;
			}
					
		// Compute the rarity value of the solution
		computeRarity(solution);
		
		// Update population
		int index = 0;
		boolean ok = false;
		SolutionWrapper dropped = solution;
		synchronized (this)
		{
			while(!ok)
			{
				Partition partition = partitions.get(index);
				boolean meetsStd = true;
				double rarity = dropped.getRarity();
				if(Double.isNaN(rarity) || rarity >= partition.getRarityStd())
				{
					dropped = partition.addSolution(dropped);
					index++;
				}
				else
					meetsStd = false;
				ok = !meetsStd || dropped == null || index == partitions.size();
			}
		}
		if(dropped != null)
		{
			if(dropped.equals(solution))
				refuseCount++;
			else
			{
				addedSinceUpdate++;
				refuseCount = 0;
			}
		}
		else
		{
			if(solution.getIcIndex() == -1)
			{
				addedSinceUpdate++;
				refuseCount = 0;
			}
		}
	}

	private void updateVarDist() 
	{		
		// Create counting structures
		ArrayList<ContSeries> discSeries = new ArrayList<ContSeries>();
		if(discVars != null)
		{
			for(int i = 0 ; i < discVars.size() ; i++)
			{
				ContSeries ser = new ContSeries();
				discSeries.add(ser);
			}
		}
		ArrayList<ContSeries> contSeries = new ArrayList<ContSeries>();
		if(contVars != null)
		{
			for(int i = 0 ; i < contVars.size() ; i++)
			{
				ContSeries ser = new ContSeries();
				contSeries.add(ser);
			}
		}
		
		// Populate series
		synchronized (this)
		{
			Partition part = getPartition(0);
			for(SolutionWrapper sol : part.getSolutions())
			{
				if(discVars != null)
				{
					ArrayList<Integer> discVals = sol.getSolution().getDiscValues();
					for(int i = 0 ; i < discVars.size() ; i++)
						discSeries.get(i).addValue(discVals.get(i));
				}
				if(contVars != null)
				{
					ArrayList<Double> contVals = sol.getSolution().getContValues();
					for(int i = 0 ; i < contVars.size() ; i++)
						contSeries.get(i).addValue(contVals.get(i));
				}
			}
		}
		
		// Adjust distributions
		if(discVars != null)
			for(int i = 0 ; i < discVars.size() ; i++)
				discVars.get(i).adjustDistribution(discSeries.get(i));
		if(contVars != null)
			for(int i = 0 ; i < contVars.size() ; i++)
				contVars.get(i).adjustDistribution(contSeries.get(i));
		
		// Re-compute rarity for solutions in population
		synchronized (this)
		{
			for(Partition part :getPartitions())
				for(SolutionWrapper sol : part.getSolutions())
					computeRarity(sol);
		}
	}
	
	@Override
	public double computeFirstPartStDev()
	{
		synchronized (this) 
		{
			ContSeries stDevs = new ContSeries();
			if (discVars != null) 
			{
				for (int i = 0; i < discVars.size(); i++) 
				{
					ContSeries values = new ContSeries();
					for (SolutionWrapper sol : getPartition(0).getSolutions())
						values.addValue(sol.getSolution().getDiscValues().get(i));
					stDevs.addValue(values.getStDev());
				}
			}
			if (contVars != null) 
			{
				for (int i = 0; i < contVars.size(); i++)
				{
					ContSeries values = new ContSeries();
					for (SolutionWrapper sol : getPartition(0).getSolutions())
						values.addValue(sol.getSolution().getContValues().get(i));
					stDevs.addValue(values.getStDev());
				}
			}
			return stDevs.getMean();
		}
	}

	/**
	 * Computes and assigns the rarity of a solution given the distributions of the values of each
	 * variable in the solutions in the first partition of the population. The rarity measures how 
	 * different or rare is the value. The rarity is a value from 0 (average solution) to 1 (rarest 
	 * solution).
	 * @param solution The solution to compute and assign its rarity
	 */
	private void computeRarity(SolutionWrapper solution)
	{
		ContSeries indRarities = new ContSeries();
		if(discVars != null)
			for(int i = 0 ; i < discVars.size() ; i++)
			{
				int value = solution.getSolution().getDiscValues().get(i);
				indRarities.addValue(discVars.get(i).computeRarity(value));
			}
		if(contVars != null)
			for(int i = 0 ; i < contVars.size() ; i++)
			{
				double value = solution.getSolution().getContValues().get(i);
				indRarities.addValue(contVars.get(i).computeRarity(value));
			}
		solution.setRarity(indRarities.getMean());
	}
	
	@Override
	public ArrayList<Solution> getAllSolutionsBestToWorst() 
	{
		ArrayList<Solution> list = new ArrayList<>();
		for(Partition part : partitions)
		{
			Iterator<SolutionWrapper> iter = part.getSolutions().descendingIterator();
			while(iter.hasNext())
				list.add(iter.next().getSolution());
		}
		return list;
	}
	
	@Override
	public void clear() 
	{
		// TODO Auto-generated method stub
	}
	
	@Override
	public boolean contains(SolutionWrapper solution) 
	{
		// TODO Auto-generated method stub
		return false;
	}

	@Override
	public int getRefuseCount() 
	{
		return refuseCount;
	}

	@Override
	public int getAddedSinceUpdate() 
	{
		return addedSinceUpdate;
	}
	
	@Override
	public ArrayList<SolutionWrapper> getAllSolutions() 
	{
		// TODO Auto-generated method stub
		return null;
	}

	@Override
	public ArrayList<SolutionWrapper> select(int count) 
	{
		// TODO Auto-generated method stub
		return null;
	}

	@Override
	public ArrayList<SolutionWrapper> select(int count, double greed) 
	{
		// TODO Auto-generated method stub
		return null;
	}

	@Override
	public HashSet<SolutionWrapper> selectSet(int count) 
	{
		// TODO Auto-generated method stub
		return null;
	}

	@Override
	public HashSet<SolutionWrapper> selectSet(int count, double greed) 
	{
		// TODO Auto-generated method stub
		return null;
	}

}
