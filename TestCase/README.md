## Repository Structure

The repository contains three OpenFOAM test cases used to reproduce the main results presented in the manuscript.

```text
TestCase/
├── Without-airInjection/
├── With-airInjection_FROM_Tap1_Q1LperMin/
└── With-airInjection_FROM_Tap5_Q1LperMin/
Test Cases and Corresponding Manuscript Results
Manuscript result	Case folder	Description
Fig. 2, Fig. 3, Table 3, Table 5	TestCase/Without-airInjection/	Baseline cavitating hydrofoil case without air injection, with Q = 0 L/min.
Fig. 4–7, Fig. 9–13, Fig. 15–17, Table 6–9	TestCase/With-airInjection_FROM_Tap1_Q1LperMin/	Main air-injection case with Q = 1 L/min, using injection through Tap 1.
Fig. 14	TestCase/With-airInjection_FROM_Tap5_Q1LperMin/	Additional validation case with Q = 1 L/min, using injection through Tap 5 for comparison with the corresponding experimental cavitation-area data.

Most air-injection results presented in the manuscript are obtained using injection through Tap 1. The Tap 5 case is provided separately because the cavitation-area validation shown in Fig. 14 follows the corresponding experimental configuration available in the reference data.
