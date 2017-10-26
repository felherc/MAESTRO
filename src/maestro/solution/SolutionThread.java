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

import java.util.ArrayList;

import maestro.MAESTROptimizer;
import maestro.solution.Solution;
import maestro.solution.SolutionRoot;
import maestro.solution.SolutionWrapper;

/**
 * This class invokes the generation and processing of new solutions. The processing of a solution
 * root involves the creation of a <code>Solution</code> object as defined by the user, and thus 
 * making it conforming to existing restrictions and computing its fitness values to be compared
 * to other solutions. The solution is wrapped within a <code>SolutionWrapper</code> object that 
 * keeps the identifier of the solution generator and its rarity value. The whole process is 
 * executed in separate threads to take advantage of of multi-core CPUs and thus speeding up run 
 * times.
 * @author Felipe Hernández
 */
public class SolutionThread extends Thread
{

	// --------------------------------------------------------------------------------------------
	// Attributes
	// --------------------------------------------------------------------------------------------
	
	/**
	 * Instance of the MAESTRO manager class
	 */
	private MAESTROptimizer optimizer;
	
	/**
	 * The initial solution object provided by the user from where to create new solutions
	 */
	private Solution solution;
	
	/**
	 * True if the thread should continue processing new solutions
	 */
	private boolean run;
	
	// --------------------------------------------------------------------------------------------
	// Constructor
	// --------------------------------------------------------------------------------------------
	
	/**
	 * Creates a new solution thread
	 * @param manager Instance of the MAESTRO manager class
	 * @param solution The initial solution object provided by the user from where to create new 
	 * solutions
	 */
	public SolutionThread(MAESTROptimizer manager, Solution solution)
	{
		this.optimizer = manager;
		this.solution = solution;
		run = true;
	}
	
	// --------------------------------------------------------------------------------------------
	// Methods
	// --------------------------------------------------------------------------------------------
	
	/**
	 * Tells the solution thread to stop processing new solutions
	 */
	public void end()
	{
		run = false;
	}
	
	/**
	 * Invokes the generation and processing of new solutions within an inner cycle run. The thread
	 * stops processing new solutions when the inner cycle ends.
	 */
	public void run()
	{
		run = true;
		while(run)
		{
			// Obtain solution root
			SolutionRoot root	= null;
			try
			{
				root			= optimizer.getSolutionRoot();
			} catch (Exception e)
			{
				e.printStackTrace();
			}
			
			// Analyze solution root
			try
			{
				if(root != null)
				{
					int icIndex = root.getIcIndex();
					int id = root.getIndex();
					int genIndex = root.getGenIndex();
					ArrayList<Integer> discValues = root.getDiscValues();
					ArrayList<Double> contValues = root.getContValues();
					Solution sol = solution.createNew(id, discValues, contValues);
					SolutionWrapper wrap = new SolutionWrapper(sol, genIndex);
					wrap.setIcIndex(icIndex);
					optimizer.addSolution(wrap);
				}
				else
					run = false;
			} catch (Exception e)
			{
				e.printStackTrace();
			}
		}
	}
	
}
