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
import maestro.gen.Generator;
import maestro.solution.Solution;
import maestro.solution.SolutionRoot;

/**
 * This wrapper class contains the objects that implement the <code>Generator</code> interface. The 
 * wrapper keeps additional information which includes the suffix of the generator to differentiate
 * it from other generators of the same type, the number of solutions generated that made it into 
 * the current population and the total number of solutions generated.
 * @author Felipe Hernández
 */
public class GenWrapper 
{

	// --------------------------------------------------------------------------------------------
	// Attributes
	// --------------------------------------------------------------------------------------------
	
	/**
	 * Suffix for the identifier of the generator. Used when more than one of the generators has 
	 * the same identifier.
	 */
	private String suffix;
	
	/**
	 * The generator being contained
	 */
	private Generator generator;
	
	/**
	 * The total number of solutions generated
	 */
	private int genTotal;
	
	// --------------------------------------------------------------------------------------------
	// Constructor
	// --------------------------------------------------------------------------------------------
	
	/**
	 * Initializes the generator wrapper
	 * @param generator The generator being contained
	 */
	public GenWrapper(Generator generator)
	{
		suffix = "";
		this.generator = generator;
		genTotal = 0;
	}
	
	// --------------------------------------------------------------------------------------------
	// Methods
	// --------------------------------------------------------------------------------------------
	
	/**
	 * @return The identifier of the generator method
	 */
	public String getId()
	{
		return generator.getId() + suffix;
	}
	
	/**
	 * @return The short identifier of the generator method
	 */
	public String getShortId()
	{
		return generator.getShortId() + suffix;
	}
	
	/**
	 * @return Suffix for the identifier of the generator. Used when more than one of the 
	 * generators has the same identifier.
	 */
	public String getSuffix()
	{
		return suffix;
	}
	
	/**
	 * @param suffix Suffix for the identifier of the generator. Used when more than one of the 
	 * generators has the same identifier.
	 */
	public void setSuffix(String suffix)
	{
		this.suffix = suffix;
	}
	
	/**
	 * @return The total number of solutions generated
	 */
	public int getGenTotal()
	{
		return genTotal;
	}
	
	/**
	 * @return A line that shows the values of the different parameters of the generator method
	 */
	public String getParamSummary()
	{
		return generator.getParamSummary();
	}
	
	/**
	 * Adds a new discrete variable to generate values for
	 * @param variable The discrete variable to add
	 */
	public void addDiscVariable(DiscVar variable)
	{
		generator.addDiscVariable(variable);
	}
	
	/**
	 * Adds a new continuous variable to generate values for
	 * @param variable The continuous variable to add
	 */
	public void addContVariable(ContVar variable)
	{
		generator.addContVariable(variable);
	}
	
	/**
	 * Deletes all the discrete and continuous variables
	 */
	public void clearVariables()
	{
		generator.clearVariables();
	}
	
	/**
	 * Generates an specific number of candidate solutions based on the current solution population
	 * @param population The current list of solutions
	 * @param number The number of new candidate solutions to be created
	 * @return A list with an specific number of generated candidate solutions
	 */
	public ArrayList<SolutionRoot> generateSolutions(ArrayList<Solution> population, int number)
	{
		ArrayList<SolutionRoot> roots = generator.generateSolutions(population, number);
		genTotal += roots.size();
		return roots;
	}
	
}
