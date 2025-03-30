package main

import "strings"

func SAISEntryPoint(s []int, K int) []int {
	n := len(s)
	return SAIS(s, K, n, make([]int, n), make([]int, n))
}

// SAIS constructs the suffix array for s using the SAIS algorithm.
// s is expected to have a trailing sentinel (0).
func SAIS(s []int, K int, n int, SA []int, lmsNames []int) []int {
	SA = SA[:n]
	for i := range SA {
		SA[i] = -1
	}
	if n == 0 {
		return SA
	}
	if n == 1 {
		SA[0] = 0
		return SA
	}
	t := make([]bool, n)
	t[n-1] = true
	for i := n - 2; i >= 0; i-- {
		if s[i] < s[i+1] {
			t[i] = true
		} else if s[i] > s[i+1] {
			t[i] = false
		} else {
			t[i] = t[i+1]
		}
	}
	var lmsPositions []int
	for i := 1; i < n; i++ {
		if t[i] && !t[i-1] {
			lmsPositions = append(lmsPositions, i)
		}
	}
	SA = induceSort(s, SA, t, K, lmsPositions)
	var sortedLMS []int
	for _, pos := range SA {
		if pos > 0 && t[pos] && !t[pos-1] {
			sortedLMS = append(sortedLMS, pos)
		}
	}
	lmsNames = lmsNames[:n]
	for i := range lmsNames {
		lmsNames[i] = -1
	}
	name := 0
	prev := -1
	for _, pos := range sortedLMS {
		if prev == -1 {
			lmsNames[pos] = name
		} else {
			if !lmsSubstringEqual(s, t, prev, pos) {
				name++
			}
			lmsNames[pos] = name
		}
		prev = pos
	}
	numNames := name + 1
	reduced := make([]int, 0, len(lmsPositions))
	for _, pos := range lmsPositions {
		reduced = append(reduced, lmsNames[pos])
	}
	var reducedSA []int
	if numNames < len(reduced) {
		reducedSA = SAIS(reduced, numNames, len(reduced), SA, lmsNames)
	} else {
		reducedSA = make([]int, len(reduced))
		for i, name := range reduced {
			reducedSA[name] = i
		}
	}
	orderedLMS := make([]int, len(reducedSA))
	for i, idx := range reducedSA {
		orderedLMS[i] = lmsPositions[idx]
	}
	for i := range SA {
		SA[i] = -1
	}
	SA = induceSort(s, SA, t, K, orderedLMS)
	return SA
}

func induceSort(s []int, SA []int, t []bool, K int, lms []int) []int {
	bs := computeBucketSizes(s, K)
	bucketTails := computeBucketTails(bs)
	for i := len(lms) - 1; i >= 0; i-- {
		pos := lms[i]
		c := s[pos]
		SA[bucketTails[c]] = pos
		bucketTails[c]--
	}
	bucketHeads := computeBucketHeads(bs)
	for i := range SA {
		pos := SA[i]
		if pos > 0 && !t[pos-1] {
			c := s[pos-1]
			SA[bucketHeads[c]] = pos - 1
			bucketHeads[c]++
		}
	}
	bucketTails = computeBucketTails(bs)
	for i := len(SA) - 1; i >= 0; i-- {
		pos := SA[i]
		if pos > 0 && t[pos-1] {
			c := s[pos-1]
			SA[bucketTails[c]] = pos - 1
			bucketTails[c]--
		}
	}
	return SA
}

func computeBucketSizes(s []int, K int) []int {
	bs := make([]int, K)
	for i := 0; i < len(s); i++ {
		bs[s[i]]++
	}
	return bs
}

func computeBucketHeads(bs []int) []int {
	heads := make([]int, len(bs))
	sum := 0
	for i, v := range bs {
		heads[i] = sum
		sum += v
	}
	return heads
}

func computeBucketTails(bs []int) []int {
	tails := make([]int, len(bs))
	sum := 0
	for i, v := range bs {
		sum += v
		tails[i] = sum - 1
	}
	return tails
}

func lmsSubstringEqual(s []int, t []bool, i, j int) bool {
	n := len(s)
	for {
		if s[i] != s[j] {
			return false
		}
		iIsLMS := (i > 0 && t[i] && !t[i-1])
		jIsLMS := (j > 0 && t[j] && !t[j-1])
		if iIsLMS && jIsLMS {
			return true
		}
		if iIsLMS != jIsLMS {
			return false
		}
		i++
		j++
		if i >= n || j >= n {
			break
		}
	}
	return false
}

// searchSequence uses binary search on the suffix entries to locate all occurrences of query.
func searchSequence(genome string, entries []SuffixEntry, query string) []SuffixEntry {
	lb := lowerBound(genome, entries, query)
	if lb == -1 {
		return []SuffixEntry{}
	}
	ub := upperBound(genome, entries, query)
	return entries[lb:ub]
}

func lowerBound(genome string, entries []SuffixEntry, query string) int {
	lo := 0
	hi := len(entries)
	for lo < hi {
		mid := (lo + hi) / 2
		if compareSuffix(genome, entries[mid].Pos, query) < 0 {
			lo = mid + 1
		} else {
			hi = mid
		}
	}
	if lo < len(entries) && strings.HasPrefix(genome[entries[lo].Pos:], query) {
		return lo
	}
	return -1
}

func upperBound(genome string, entries []SuffixEntry, query string) int {
	lo := 0
	hi := len(entries)
	for lo < hi {
		mid := (lo + hi) / 2
		if compareSuffix(genome, entries[mid].Pos, query) <= 0 {
			lo = mid + 1
		} else {
			hi = mid
		}
	}
	return lo
}

func compareSuffix(genome string, pos int, query string) int {
	i := 0
	for i < len(query) && pos+i < len(genome) {
		if genome[pos+i] != query[i] {
			return int(genome[pos+i]) - int(query[i])
		}
		i++
	}
	if i < len(query) {
		return -1
	}
	return 0
}

// encodeString converts the input string into an int slice and appends a sentinel (0).
// Each character is shifted by +1 so that 0 is reserved as the sentinel.
func encodeString(s string) ([]int, int) {
	n := len(s)
	encoded := make([]int, n+1)
	maxVal := 0
	for i, ch := range s {
		encoded[i] = int(ch) + 1
		if encoded[i] > maxVal {
			maxVal = encoded[i]
		}
	}
	encoded[n] = 0
	return encoded, maxVal + 1
}
