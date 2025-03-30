package main

import (
	"reflect"
	"strings"
	"testing"
)

func TestSuffixArray(t *testing.T) {
	// testCases holds the input string and its expected suffix array.
	// The expected suffix array is based on the encoded string that includes the trailing sentinel.
	// For example:
	// For "banana", the expected suffix array is [6, 5, 3, 1, 0, 4, 2] (index 6 corresponds to the sentinel).
	testCases := []struct {
		input    string
		expected []int
	}{
		{
			input:    "banana",
			expected: []int{6, 5, 3, 1, 0, 4, 2},
		},
		{
			input:    "mississippi",
			expected: []int{11, 10, 7, 4, 1, 0, 9, 8, 6, 3, 5, 2},
		},
		{
			input:    "a",
			expected: []int{1, 0},
		},
		{
			input:    "",
			expected: []int{0},
		},
	}

	for _, tc := range testCases {
		encoded, alphabetSize := encodeString(tc.input)
		sa := SAISEntryPoint(encoded, alphabetSize)
		if !reflect.DeepEqual(sa, tc.expected) {
			t.Errorf("Suffix array for %q: got %v, expected %v", tc.input, sa, tc.expected)
		}
	}
}

func BenchmarkSAIS(b *testing.B) {
	// Generate a synthetic genome sequence.
	// This simulates a relatively large input for benchmarking.
	genome := strings.Repeat("ACGT", 10000000) // 40,000,000 characters
	encoded, alphabetSize := encodeString(genome)

	// Reset the timer so that the setup time (string generation & encoding) is not included.
	b.ResetTimer()

	for i := 0; i < b.N; i++ {
		// Run the SAIS algorithm for each iteration.
		_ = SAISEntryPoint(encoded, alphabetSize)
	}
}
