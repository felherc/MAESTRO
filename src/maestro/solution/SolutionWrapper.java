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

package maestro.solution;

import maestro.pop.Partition;
import maestro.solution.Solution;
import maestro.solution.SolutionWrapper;

/**
 * This wrapper class contains the objects that implement the <code>Solution</code> interface. The 
 * wrapper keeps additional information which includes the index of the inner cycle run in which 
 * the solution was created, the partition the solution belongs to, the identifier of the algorithm 
 * that generated the solution and the rarity of the solution. Additionally, this class implements 
 * a <code>Comparable</code> interface that allows solutions to be sorted within the population 
 * partition using both the actual fitness function and the rarity standards of the population 
 * partitions: solutions that meet the rarity standard are preferred over those that don't, and two 
 * solutions that meet the standard are compared through their fitness values.
 * @author Felipe Hernández
 */
public class SolutionWrapper implements Comparable<SolutionWrapper>
{

	// --------------------------------------------------------------------------------------------
	// Attributes
	// --------------------------------------------------------------------------------------------
	
	/**
	 * The index of the inner cycle run in which the solution was created
	 */
	private int icIndex;
	
	/**
	 * The population partition that allocates the solution
	 */
	private Partition partition;
	
	/**
	 * The solution being wrapped
	 */
	private Solution solution;
	
	/**
	 * The index of the generator algorithm that created the solution
	 */
	private int genIndex;
	
	/**
	 * The measure of how different or rare is the solution in comparison with the solutions in the 
	 * first partition of the population. The rarity is a value form 0 (average solution) to 1 
	 * (rarest solution) and is computed as the average rarity of the solution within each of the 
	 * variables of the problem. Depending of the type of variable, the rarity is computed 
	 * differently.
	 */
	private double rarity;
	
	// --------------------------------------------------------------------------------------------
	// Constructors
	// --------------------------------------------------------------------------------------------
	
	/**
	 * Initializes the solution wrapper object
	 * @param solution The solution to be wrapped
	 * @param genIndex The index of the generator algorithm that created the solution
	 * @param rarity The measure of how different or rare is the solution in comparison with all 
	 * the solutions in the current population
	 */
	public SolutionWrapper(Solution solution, int genIndex, double rarity)
	{
		partition = null;
		this.solution = solution;
		this.genIndex = genIndex;
		this.rarity = rarity;
	}
	
	/**
	 * Initializes the solution wrapper object
	 * @param solution The solution to be wrapped
	 * @param genIndex The index of the generator algorithm that created the solution
	 */
	public SolutionWrapper(Solution solution, int genIndex)
	{
		partition = null;
		this.solution = solution;
		this.genIndex = genIndex;
		this.rarity = Double.NaN;
	}
	
	// --------------------------------------------------------------------------------------------
	// Methods
	// --------------------------------------------------------------------------------------------
	
	/**
	 * @return The identifier of the solution being wrapped
	 */
	public String getId()
	{
		return solution.getId();
	}
	
	/**
	 * @return The index of the inner cycle run in which the solution was created
	 */
	public int getIcIndex() 
	{
		return icIndex;
	}

	/**
	 * @param icIndex The index of the inner cycle run in which the solution was created
	 */
	public void setIcIndex(int icIndex) 
	{
		this.icIndex = icIndex;
	}
	
	/**
	 * @return The population partition that allocates the solution
	 */
	public Partition getPartition() 
	{
		return partition;
	}

	/**
	 * @param partition The population partition that allocates the solution
	 */
	public void setPartition(Partition partition) 
	{
		this.partition = partition;
	}

	/**
	 * @return The solution being wrapped
	 */
	public Solution getSolution() 
	{
		return solution;
	}

	/**
	 * @return The index of the generator algorithm that created the solution
	 */
	public int getGenIndex() 
	{
		return genIndex;
	}
	
	/**
	 * @return The measure of how different or rare is the solution in comparison with all the
	 * solutions in the current population
	 */
	public double getRarity() 
	{
		return rarity;
	}

	/**
	 * @param rarity The measure of how different or rare is the solution in comparison with all the
	 * solutions in the current population
	 */
	public void setRarity(double rarity) 
	{
		this.rarity = rarity;
	}
	
	@Override
	public int compareTo(SolutionWrapper other) 
	{
		if(solution.getId().equals(other.getId()))
			return 0;
		
		double rarityStd = partition.getRarityStd();
		double otherRarity = other.getRarity();
		if(!Double.isNaN(rarity) && !Double.isNaN(otherRarity))
		{
			boolean okRarity = rarity >= rarityStd;
			boolean okRarityOther = otherRarity >= rarityStd;
			if(okRarity)
				if(okRarityOther)
					return solution.compareTo(other.getSolution());
				else
					return 1;
			else
				if(okRarityOther)
					return -1;
				else
					if(rarity == otherRarity)
						return solution.compareTo(other.getSolution());
					else
						return rarity > otherRarity ? 1 : -1;
		}
		else
			return solution.compareTo(other.getSolution());
	}
	
}