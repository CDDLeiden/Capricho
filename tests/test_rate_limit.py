import time
import unittest
from concurrent.futures import ThreadPoolExecutor, as_completed
from contextlib import redirect_stderr
from io import StringIO

from Capricho.core.rate_limit import rate_limit
from Capricho.logger import setup_logger


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
            if "starting execution at" in line:
                time_str = line.split("starting execution at")[1].split(".")[0]
                execution_times.append(float(time_str))

        self.assertGreaterEqual(len(execution_times), num_calls, "Not enough execution start times logged")
        # Sort execution times if necessary, though as_completed and single lock should order them
        execution_times.sort()

        # Calculate time differences between executions
        # Only relevant for (num_calls - 1) intervals
        time_diffs = [execution_times[i + 1] - execution_times[i] for i in range(len(execution_times) - 1)]
        avg_time_diff = sum(time_diffs) / len(time_diffs) if time_diffs else 0

        # Assertions
        self.assertEqual(len(results), num_calls, "All function calls should have completed")

        # The total time should be approximately (num_calls - 1) * min_interval + (time for last call)
        # For 40 calls and 0.2s interval, this is 39 * 0.2 = 7.8s plus function execution time.
        # The current assertion (num_calls / 5) - 0.2 is (40/5) - 0.2 = 8 - 0.2 = 7.8s, which is correct.
        self.assertGreaterEqual(
            total_time,
            (num_calls / 5) - 0.2,
            "Total time should be at least (num_calls / max_per_second) minus one interval",
        )
        self.assertAlmostEqual(
            avg_time_diff,
            0.2,
            delta=0.02,  # A slightly larger delta might still be reasonable given OS scheduling variations
            msg="Average time between executions should be close to 0.2 seconds",
        )

        # Check for rate limit messages in log
        rate_limit_messages = [line for line in log_output.split("\n") if "Rate limit exceeded" in line]
        self.assertGreaterEqual(
            len(rate_limit_messages),
            num_calls - 1,
            "There should be appropriate rate limit messages in the log",
        )

        print(
            f"Test completed. Total time: {total_time:.2f}s, Average time between executions: {avg_time_diff:.4f}s"
        )


if __name__ == "__main__":
    unittest.main()
