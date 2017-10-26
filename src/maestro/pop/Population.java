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

import maestro.solution.Solution;
import maestro.solution.SolutionWrapper;

/**
 * Interfaces implementations of partitioned solution populations for single-objective optimization 
 * problems.
 * @author Felipe Hernández
 */
public interface Population
{

	// --------------------------------------------------------------------------------------------
	// Methods
	// --------------------------------------------------------------------------------------------
	
	/**
	 * Clears the population of solutions
	 */
	public void clear();
	
	/**
	 * @returns The number of solutions currently in the population
	 */
	public int size();
	
	/**
	 * @return The target size of the current population
	 */
	public int getCapacity();
	
	/**
	 * @param capacity The target size of the current population
	 */
	public void setCapacity(int capacity);
	
	/**
	 * Returns a partition given the index. Partitions are indexed from 0 to the number of 
	 * partitions minus one.
	 * @param index
	 */
	public Partition getPartition(int index);
	
	/**
	 * @return The list of all the partitions in the population
	 */
	public ArrayList<Partition> getPartitions();
	
	/**
	 * @return The number of new solutions refused from the population since the last solution was 
	 * accepted
	 */
	public int getRefuseCount();
	
	/**
	 * @return The number of processed solutions added to the population since the last update of 
	 * the distributions of the variables
	 */
	public int getAddedSinceUpdate();
	
	/**
	 * @return The percentage of new solutions accepted compared to the size of the population 
	 * before the distributions of the variables are updated
	 */
	public double getUpdateDist();
	
	/**
	 * @param updateDist The percentage of new solutions accepted compared to the size of the 
	 * population before the distributions of the variables are updated
	 */
	public void setUpdateDist(double updateDist);
	
	/**
	 * @param solution The solution to analyze
	 * @return True if an identical solution already exists in the population. That is, if the 
	 * defining values of the provided solution are the same that those in any solution in the
	 * population.
	 */
	public boolean contains(SolutionWrapper solution);
	
	/**
	 * Forcefully adds a solution to the population without taking into account its fitness and
	 * rarity
	 * @param solution The solution to add
	 */
	public void forceAddSolution(SolutionWrapper solution);
	
	/**
	 * Offers a new solution to be considered for addition to the population
	 * @param solution the solution offered
	 */
	public void offerSolution(SolutionWrapper solution);
	
	/**
	 * @return A list with all the solutions in the population
	 */
	public ArrayList<SolutionWrapper> getAllSolutions();
	
	/**
	 * @return A list with all the solutions in the population arranging them from the first 
	 * partition to the last, and from the most fitted to the least fitted in each partition.
	 * Called to produce new solution roots.
	 */
	public ArrayList<Solution> getAllSolutionsBestToWorst();
	
	/**
	 * @param count The number of requested solutions
	 * @return A list of solutions in the population. A custom strategy can be implemented.
	 */
	public ArrayList<SolutionWrapper> select(int count);
	
	/**
	 * @param count The number of requested solutions
	 * @param greed A value between -1.0 and 1.0 that represents the expected mean "quality" of the 
	 * returned solutions. 1.0 if the selected solutions should be only among the best according to 
	 * the objectives. -1.0 if they should be only among the worst. 0.0 if there is no preference.
	 * @return A list of solutions in the population. A custom strategy can be implemented.
	 */
	public ArrayList<SolutionWrapper> select(int count, double greed);
	
	/**
	 * @param count The number of requested solutions
	 * @return A set of solutions in the population. The number of solutions returned may be 
	 * smaller than the number requested. A custom strategy can be implemented.
	 */
	public HashSet<SolutionWrapper> selectSet(int count);
	
	/**
	 * @param count The number of requested solutions
	 * @param greed A value between -1.0 and 1.0 that represents the expected mean "quality" of the 
	 * returned solutions. 1.0 if the selected solutions should be only among the best according to 
	 * the objectives. -1.0 if they should be only among the worst. 0.0 if there is no preference.
	 * @return A set of solutions in the population. The number of solutions returned may be 
	 * smaller than the number requested. A custom strategy can be implemented.
	 */
	public HashSet<SolutionWrapper> selectSet(int count, double greed);
	
	/**
	 * @return The average standard deviation of the values of the solutions in the first 
	 * partition of the population
	 */
	public double computeFirstPartStDev();
	
}
