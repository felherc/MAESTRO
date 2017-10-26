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

package maestro;

/**
 * This class stores the information of an inner cycle run
 * @author Felipe Hernández
 */
public class ICLogEntry 
{
	
	// --------------------------------------------------------------------------------------------  
	// Constants
	// --------------------------------------------------------------------------------------------
	
	/**
	 * Identifier of the termination criterion if the maximum number of refused solutions was 
	 * exceeded
	 */
	public final static int REFUSED = 1;
	
	/**
	 * Identifier of the termination criterion if the uniformity of the first partition of the 
	 * population was less than the minimum
	 */
	public final static int UNIFORMITY = 2;
	
	/**
	 * Identifier of the termination criterion if the whole optimization process ended either 
	 * because the time limit or the maximum number of solutions was reached
	 */
	public final static int OC_TERMINATED = 3;

	// --------------------------------------------------------------------------------------------  
	// Attributes
	// --------------------------------------------------------------------------------------------

	/**
	 * The size of the population for the cycle
	 */
	private int population;
	
	/**
	 * The number of solutions processed within the cycle
	 */
	private int solutions;
	
	/**
	 * The accumulated number of solutions processed up to this cycle
	 */
	private int accumSolutions;
	
	/**
	 * True if a better solution than the best-so-far was found during the cycle. False by default.
	 */
	private boolean foundBest;
	
	/**
	 * The reason why the cycle was terminated as defined by the class constants:<ul>
	 * <li><code>REFUSED</code>: the maximum number of refused solutions was exceeded
	 * <li><code>UNIFORMITY</code>: the uniformity of the first partition of the population was 
	 * less than the minimum
	 * <li><code>OC_TERMINATED</code>: the whole optimization process ended either because the time 
	 * limit or the maximum number of solutions was reached</ul>
	 */
	private int termination;
	
	// --------------------------------------------------------------------------------------------  
	// Constructor
	// --------------------------------------------------------------------------------------------
	
	/**
	 * Creates a new log entry
	 * @param population The size of the population for the cycle
	 */
	public ICLogEntry(int population)
	{
		this.population = population;
		solutions = 0;
		accumSolutions = 0;
		foundBest = false;
		termination = OC_TERMINATED;
	}
	
	// --------------------------------------------------------------------------------------------  
	// Methods
	// --------------------------------------------------------------------------------------------
	
	/**
	 * @return The size of the population for the cycle
	 */
	public int getPopulation() 
	{
		return population;
	}

	/**
	 * @param population The size of the population for the cycle
	 */
	public void setPopulation(int population) 
	{
		this.population = population;
	}
	
	/**
	 * @return The number of solutions processed within the cycle
	 */
	public int getSolutions() 
	{
		return solutions;
	}

	/**
	 * @param solutions The number of solutions processed within the cycle
	 */
	public void setSolutions(int solutions) 
	{
		this.solutions = solutions;
	}
	
	
	/**
	 * @return The accumulated number of solutions processed up to this cycle
	 */
	public int getAccumSolutions() 
	{
		return accumSolutions;
	}

	/**
	 * @param accumSolutions The accumulated number of solutions processed up to this cycle
	 */
	public void setAccumSolutions(int accumSolutions) 
	{
		this.accumSolutions = accumSolutions;
	}

	/**
	 * @return True if a better solution than the best-so-far was found during the cycle
	 */
	public boolean foundBest() 
	{
		return foundBest;
	}

	/**
	 * @param foundBest True if a better solution than the best-so-far was found during the cycle
	 */
	public void setFoundBest(boolean foundBest) 
	{
		this.foundBest = foundBest;
	}

	/**
	 * @return The reason why the cycle was terminated as defined by the class constants:<ul>
	 * <li><code>REFUSED</code>: the maximum number of refused solutions was exceeded
	 * <li><code>UNIFORMITY</code>: the uniformity of the first partition of the population was 
	 * less than the minimum
	 * <li><code>OC_TERMINATED</code>: the whole optimization process ended either because the time 
	 * limit or the maximum number of solutions was reached</ul>
	 */
	public int getTermination() 
	{
		return termination;
	}

	/**
	 * @param termination The reason why the cycle was terminated as defined by the class 
	 * constants:<ul>
	 * <li><code>REFUSED</code>: the maximum number of refused solutions was exceeded
	 * <li><code>UNIFORMITY</code>: the uniformity of the first partition of the population was 
	 * less than the minimum
	 * <li><code>OC_TERMINATED</code>: the whole optimization process ended either because the time 
	 * limit or the maximum number of solutions was reached</ul>
	 */
	public void setTermination(int termination) 
	{
		this.termination = termination;
	}
	
}
