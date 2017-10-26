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
 * This interface allows MAESTRO to communicate with the caller application
 * @author Felipe Hernández
 */
public interface Monitor 
{

	// --------------------------------------------------------------------------------------------
	// Methods
	// --------------------------------------------------------------------------------------------
	
	/**
	 * Notifies of the termination of the optimization process, allowing the calling of the 
	 * reporting methods in the <code>MAESTROptimizer</code> class.
	 * @param message The message with the reason for termination
	 */
	public void terminate(String message);
	
	/**
	 * Resets the local variables to allow the beginning of the execution of a new optimization 
	 * process. This method is called when the <code>reset()</code> method in the 
	 * <code>MAESTROptimizer</code> class is called.
	 */
	public void reset();
	
}
