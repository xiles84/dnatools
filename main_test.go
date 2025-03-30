package main

import (
	"os"
	"strings"
	"testing"
)

func TestApplicationOutput(t *testing.T) {
	// Create a temporary directory for test files.
	tempDir := t.TempDir()

	// Create a temporary genome file.
	// This file contains two DNA sequences, one per line.
	genomeContent := "ACGT\nTGCA\n"
	genomeFile := tempDir + "/genome_test.txt"
	if err := os.WriteFile(genomeFile, []byte(genomeContent), 0644); err != nil {
		t.Fatalf("Failed to write genome file: %v", err)
	}

	// Create a temporary pattern file for trie search.
	// Each pattern is on its own line.
	patternContent := "ACG\nTGC\n"
	patternFile := tempDir + "/patterns_test.txt"
	if err := os.WriteFile(patternFile, []byte(patternContent), 0644); err != nil {
		t.Fatalf("Failed to write pattern file: %v", err)
	}

	// Save original os.Args and os.Stdout to restore later.
	origArgs := os.Args
	origStdout := os.Stdout

	// Set up command-line arguments for trie search mode.
	// -f for genome file and -t for pattern file.
	os.Args = []string{"cmd", "-f", genomeFile, "-t", patternFile}

	// Capture standard output.
	r, w, err := os.Pipe()
	if err != nil {
		t.Fatalf("Failed to create pipe: %v", err)
	}
	os.Stdout = w

	// Run the application.
	// (Note: main() uses flag.Parse so this should be run once per test case.)
	main()

	// Close writer and restore os.Stdout.
	w.Close()
	os.Stdout = origStdout

	// Read the captured output.
	var outputBuilder strings.Builder
	buf := make([]byte, 1024)
	for {
		n, err := r.Read(buf)
		if n > 0 {
			outputBuilder.Write(buf[:n])
		}
		if err != nil {
			break
		}
	}
	output := outputBuilder.String()

	// Check that the output contains the expected annotated positions.
	// The concatenated genome is "ACGT$TGCA", so:
	// - "ACG" should be found at position 0 on DNA line 0.
	// - "TGC" should be found at position 5 on DNA line 1.
	expectedSubstr1 := `Pattern "ACG" found at positions: [(0, line 0)]`
	expectedSubstr2 := `Pattern "TGC" found at positions: [(5, line 1)]`
	if !strings.Contains(output, expectedSubstr1) {
		t.Errorf("Output does not contain expected substring for ACG. Got:\n%s", output)
	}
	if !strings.Contains(output, expectedSubstr2) {
		t.Errorf("Output does not contain expected substring for TGC. Got:\n%s", output)
	}

	// Restore original os.Args.
	os.Args = origArgs
}

// TestIndexModeOutput creates a temporary genome file, runs the application in index mode (-m),
// and checks that the output contains the "Index built" message and that the index file (sa.idx) is created.
func TestIndexModeOutput(t *testing.T) {
	// Create a temporary directory and genome file.
	tempDir := t.TempDir()
	genomeContent := "banana\n" // single DNA sequence
	genomeFile := tempDir + "/genome_index.txt"
	if err := os.WriteFile(genomeFile, []byte(genomeContent), 0644); err != nil {
		t.Fatalf("Failed to write genome file: %v", err)
	}

	origArgs := os.Args
	origStdout := os.Stdout

	// Set command-line arguments for index mode.
	os.Args = []string{"cmd", "-m", "-f", genomeFile}

	// Capture the standard output.
	r, w, err := os.Pipe()
	if err != nil {
		t.Fatalf("Failed to create pipe: %v", err)
	}
	os.Stdout = w

	// Run the application.
	main()

	// Close the writer and restore stdout.
	w.Close()
	os.Stdout = origStdout

	// Read the captured output.
	var outputBuilder strings.Builder
	buf := make([]byte, 1024)
	for {
		n, err := r.Read(buf)
		if n > 0 {
			outputBuilder.Write(buf[:n])
		}
		if err != nil {
			break
		}
	}
	output := outputBuilder.String()

	// Check that the output contains the expected index built message.
	if !strings.Contains(output, "Index built and saved to sa.idx") {
		t.Errorf("Expected index built message, got output: %s", output)
	}

	// Verify that the index file "sa.idx" was created.
	if _, err := os.Stat("sa.idx"); os.IsNotExist(err) {
		t.Errorf("Index file sa.idx was not created")
	}
	// Clean up the index file.
	os.Remove("sa.idx")
	os.Args = origArgs
}

// TestSearchModeOutput builds an index for a simple genome file and then runs a search (-s) mode.
// It verifies that the output contains the expected search results with global positions and DNA line numbers.
func TestSearchModeOutput(t *testing.T) {
	tempDir := t.TempDir()
	// Use a known genome (a single DNA sequence) that will be processed as one line.
	genomeContent := "banana\n"
	genomeFile := tempDir + "/genome_search.txt"
	if err := os.WriteFile(genomeFile, []byte(genomeContent), 0644); err != nil {
		t.Fatalf("Failed to write genome file: %v", err)
	}

	origArgs := os.Args
	origStdout := os.Stdout

	// First, build the index.
	os.Args = []string{"cmd", "-m", "-f", genomeFile}
	_, w1, err := os.Pipe()
	if err != nil {
		t.Fatalf("Failed to create pipe: %v", err)
	}
	os.Stdout = w1
	main() // build index mode
	w1.Close()
	os.Stdout = origStdout

	// Now, run the search mode for the query "ana".
	os.Args = []string{"cmd", "-s", "ana", "-f", genomeFile}
	r2, w2, err := os.Pipe()
	if err != nil {
		t.Fatalf("Failed to create pipe: %v", err)
	}
	os.Stdout = w2
	main() // search mode
	w2.Close()
	os.Stdout = origStdout

	// Read the captured output.
	var outputBuilder strings.Builder
	buf := make([]byte, 1024)
	for {
		n, err := r2.Read(buf)
		if n > 0 {
			outputBuilder.Write(buf[:n])
		}
		if err != nil {
			break
		}
	}
	output := outputBuilder.String()

	// Expected: "Sequence found at positions:" message, and for "banana" with query "ana",
	// the substring "ana" occurs at positions 1 and 3 (with DNA line 0).
	if !strings.Contains(output, "Sequence found at positions (global position, DNA line):") {
		t.Errorf("Expected search result message, got output: %s", output)
	}
	if !strings.Contains(output, "(1, 0)") || !strings.Contains(output, "(3, 0)") {
		t.Errorf("Expected search result positions for query 'ana', got output: %s", output)
	}

	// Clean up the index file and restore os.Args.
	os.Remove("sa.idx")
	os.Args = origArgs
}

func BenchmarkSAIS(b *testing.B) {
	// Generate a synthetic genome sequence for benchmarking.
	genome := strings.Repeat("ACGT", 10000000) // 40,000,000 characters
	encoded, alphabetSize := encodeString(genome)
	b.ResetTimer()
	for i := 0; i < b.N; i++ {
		_ = SAISEntryPoint(encoded, alphabetSize)
	}
}
