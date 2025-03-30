package main

import (
	"strings"
	"testing"
)

func BenchmarkSAIS(b *testing.B) {
	// Generate a synthetic genome sequence for benchmarking.
	genome := strings.Repeat("ACGT", 10000000) // 40,000,000 characters
	encoded, alphabetSize := encodeString(genome)
	b.ResetTimer()
	for i := 0; i < b.N; i++ {
		_ = SAISEntryPoint(encoded, alphabetSize)
	}
}
