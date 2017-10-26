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

/**
 * This class contains the generated discrete and continuous values of a candidate solution to be
 * analyzed
 * @author Felipe Hernández
 */
public class SolutionRoot 
{

	// --------------------------------------------------------------------------------------------
	// Attributes
	// --------------------------------------------------------------------------------------------
	
	/**
	 * The index of the solution
	 */
	private int index;
	
	/**
	 * The index of the inner cycle run in which the solution was created
	 */
	private int icIndex;
	
	/**
	 * The index of the generator algorithm that created the solution
	 */
	private int genIndex;
	
	/**
	 * List with the values of the solution for the discrete variables
	 */
	private ArrayList<Integer> discValues;
	
	/**
	 * List with the values of the solution for the continuous variables
	 */
	private ArrayList<Double> contValues;
	
	// --------------------------------------------------------------------------------------------
	// Constructor
	// --------------------------------------------------------------------------------------------
	
	/**
	 * Creates a new solution root
	 * @param discValues List with the values of the solution for the discrete variables
	 * @param contValues List with the values of the solution for the continuous variables
	 */
	public SolutionRoot(ArrayList<Integer> discValues, ArrayList<Double> contValues)
	{
		index = -1;
		icIndex = -1;
		genIndex = -1;
		this.discValues = discValues;
		this.contValues = contValues;
	}
	
	// --------------------------------------------------------------------------------------------
	// Methods
	// --------------------------------------------------------------------------------------------
	
	/**
	 * @return The index of the solution
	 */
	public int getIndex() 
	{
		return index;
	}

	/**
	 * @param index The index of the solution
	 */
	public void setIndex(int index) 
	{
		this.index = index;
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
	 * @return The index of the generator algorithm that created the solution
	 */
	public int getGenIndex() 
	{
		return genIndex;
	}

	/**
	 * @param genIndex The index of the generator algorithm that created the solution
	 */
	public void setGenIndex(int genIndex) 
	{
		this.genIndex = genIndex;
	}

	/**
	 * @return List with the values of the solution for the discrete variables
	 */
	public ArrayList<Integer> getDiscValues() 
	{
		return discValues;
	}

	/**
	 * @return List with the values of the solution for the continuous variables
	 */
	public ArrayList<Double> getContValues() 
	{
		return contValues;
	}
	
}