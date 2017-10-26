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

package maestro.gen.gd;

import maestro.gen.gd.IndexDoubleSorter;

/**
 * This class allows ordering indexes according to a double criteria
 * @author Felipe Hernández
 */
public class IndexDoubleSorter implements Comparable<IndexDoubleSorter>
{

	// --------------------------------------------------------------------------------------------
	// Attributes
	// --------------------------------------------------------------------------------------------
	
	/**
	 * The index to be ordered
	 */
	private int index;
	
	/**
	 * The value used for sorting
	 */
	private double sortingValue;
	
	// --------------------------------------------------------------------------------------------
	// Constructor
	// --------------------------------------------------------------------------------------------
	
	/**
	 * @param index 		{@link #index}
	 * @param sortingValue	{@link #sortingValue}
	 */
	public IndexDoubleSorter(int index, double sortingValue)
	{
		this.index			= index;
		this.sortingValue	= sortingValue;
	}
	
	// --------------------------------------------------------------------------------------------
	// Methods
	// --------------------------------------------------------------------------------------------

	/**
	 * @return {@link #index}
	 */
	public int getIndex() 
	{
		return index;
	}

	/**
	 * @param index {@link #index}
	 */
	public void setIndex(int index) 
	{
		this.index = index;
	}

	/**
	 * @return {@link #sortingValue}
	 */
	public double getSortingValue() 
	{
		return sortingValue;
	}

	/**
	 * @param sortingValue {@link #sortingValue}
	 */
	public void setSortingValue(double sortingValue) 
	{
		this.sortingValue = sortingValue;
	}

	@Override
	public int compareTo(IndexDoubleSorter other) 
	{
		Double d1	= new Double(sortingValue);
		Double d2	= new Double(other.getSortingValue());
		return d1.compareTo(d2);
	}	
	
}
