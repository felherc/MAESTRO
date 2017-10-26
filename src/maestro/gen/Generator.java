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

package maestro.gen;

import java.util.ArrayList;

import maestro.ContVar;
import maestro.DiscVar;
import maestro.solution.Solution;
import maestro.solution.SolutionRoot;

/**
 * This interface allows the core MAESTRO module to communicate with the different solution 
 * generator methods available to generate values for new solutions
 * @author Felipe Hernández
 */
public interface Generator 
{

	// --------------------------------------------------------------------------------------------
	// Methods
	// --------------------------------------------------------------------------------------------
	
	/**
	 * @return The identifier of the generator method
	 */
	public String getId();
	
	/**
	 * @return The short version of the identifier of the generator method
	 */
	public String getShortId();
	
	/**
	 * @return A line that shows the values of the different parameters of the generator method
	 */
	public String getParamSummary();
	
	/**
	 * Adds a new discrete variable to generate values for
	 * @param variable The discrete variable to add
	 */
	public void addDiscVariable(DiscVar variable);
	
	/**
	 * Adds a new continuous variable to generate values for
	 * @param variable The continuous variable to add
	 */
	public void addContVariable(ContVar variable);
	
	/**
	 * Deletes all the discrete and continuous variables
	 */
	public void clearVariables();
	
	/**
	 * Generates an specific number of candidate solutions based on the current solution 
	 * population. Concurrent calls to this method may occur.
	 * @param population The current ordered list of solutions starting with the most fitted
	 * @param number The number of new candidate solutions to be created
	 * @return A list with an specific number of generated candidate solutions
	 */
	public ArrayList<SolutionRoot> generateSolutions(ArrayList<Solution> population, int number);
	
}
