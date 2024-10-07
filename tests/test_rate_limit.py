import time
import unittest
from concurrent.futures import ThreadPoolExecutor, as_completed
from contextlib import redirect_stderr
from io import StringIO

from CompoundMapper.chembl.rate_limit import rate_limit
from CompoundMapper.logger import setup_logger


class TestRateLimitDecorator(unittest.TestCase):
    def setUp(self):
        self.log_capture = StringIO()
        setup_logger(level="TRACE", out_file=self.log_capture)

    @rate_limit(max_per_second=5)
    def sample_function(self, i):
        return f"Executed {i}"

    def test_rate_limit_decorator(self):
        num_calls = 40
        start_time = time.perf_counter()

        with redirect_stderr(StringIO()):  # Redirect stderr to avoid cluttering the test output
            with ThreadPoolExecutor(max_workers=10) as executor:
                futures = [executor.submit(self.sample_function, i) for i in range(num_calls)]
                results = [future.result() for future in as_completed(futures)]

        end_time = time.perf_counter()
        total_time = end_time - start_time

        # Analyze log output
        log_output = self.log_capture.getvalue()
        execution_times = []
        for line in log_output.split("\n"):
            if "Updated last_called to" in line:
                time_str = line.split("Updated last_called to")[1].split(".")[0]
                execution_times.append(float(time_str))

        # Calculate time differences between executions
        time_diffs = [execution_times[i + 1] - execution_times[i] for i in range(len(execution_times) - 1)]
        avg_time_diff = sum(time_diffs) / len(time_diffs)

        # Assertions
        self.assertEqual(len(results), num_calls, "All function calls should have completed")
        self.assertGreaterEqual(
            total_time, (num_calls / 5) - 0.2, "Total time should be at least (num_calls / max_per_second)"
        )
        self.assertAlmostEqual(
            avg_time_diff,
            0.2,
            delta=0.02,
            msg="Average time between executions should be close to 0.2 seconds",
        )

        # Check for rate limit messages in log
        rate_limit_messages = [line for line in log_output.split("\n") if "Rate limit exceeded" in line]
        self.assertGreater(len(rate_limit_messages), 0, "There should be some rate limit messages in the log")

        print(
            f"Test completed. Total time: {total_time:.2f}s, Average time between executions: {avg_time_diff:.4f}s"
        )


if __name__ == "__main__":
    unittest.main()
