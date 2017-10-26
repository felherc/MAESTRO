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

import maestro.MAESTROptimizer;

/**
 * This thread invokes the creation of new solutions on a different thread not to interrupt the 
 * processing of solutions
 * @author Felipe Hernández
 */
public class PopulatorThread extends Thread 
{

	// --------------------------------------------------------------------------------------------
	// Attributes
	// --------------------------------------------------------------------------------------------
	
	/**
	 * Instance of the MAESTRO manager class
	 */
	private MAESTROptimizer manager;
	
	// --------------------------------------------------------------------------------------------
	// Constructor
	// --------------------------------------------------------------------------------------------
	
	/**
	 * Creates a new populator thread
	 * @param manager Instance of the MAESTRO manager class
	 */
	public PopulatorThread(MAESTROptimizer manager)
	{
		this.manager = manager;
	}
	
	// --------------------------------------------------------------------------------------------
	// Method
	// --------------------------------------------------------------------------------------------
	
	/**
	 * Invokes the creation of new solutions
	 */
	public void run()
	{
		manager.populateBuffer();
	}
	
}
