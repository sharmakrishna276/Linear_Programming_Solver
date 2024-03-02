### Running the LP Tester

To test the LP (Linear Programming) solver, follow these steps:

1. **Generate Testcase and Commands Files:**
   
   Run the following command to generate the `testcase.txt` and `commands.txt` files:

   ```bash
   python tester.py > commands.txt
   ```

   - `testcase.txt`: This file contains the generated test case in the specified format as per the assignment.
   - `commands.txt`: This file contains the processing done by the PuLP module.

2. **Run LP Solver:**
   
   After generating the necessary files, execute the LP solver script:

   ```bash
   python test_later.py
   ```

   This script will read the `commands.txt` file and generate an `output.txt` file.

3. **Review Output:**

   The `output.txt` file contains the status of the LP problem and the optimal vector and costs if the status is optimal.

### File Descriptions:

- **`tester.py`:**
  Python script responsible for generating the test case and commands files.

- **`test_later.py`:**
  Python script that reads the commands file and generates the output file.

- **`testcase.txt`:**
  File containing the generated test case in the specified format.

- **`commands.txt`:**
  File containing the processing done by the PuLP module.

- **`output.txt`:**
  File containing the status of the LP problem and the optimal vector and costs if the status is optimal.

Ensure that you have Python installed on your system and the required libraries (e.g., PuLP) installed before running these scripts.