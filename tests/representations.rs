extern crate alignment_free_methods;

#[cfg(test)]
mod tests {
    use pretty_assertions::assert_eq;
    use anyhow::Result;
    
    use std::collections::HashMap;

    use alignment_free_methods::utils::representation_methods;

    #[test]
    fn test_generate_kmers() -> Result<()> {
        let kmers = representation_methods::generate_kmers(
            &vec!['A', 'C', 'T', 'G'],
            2
        )?;
        let expected = vec![&['A', 'C'], &['C', 'T'], &['T', 'G']];
        assert_eq!(kmers, expected);
        Ok(())
    }

    #[test]
    fn test_generate_kmers_with_step() -> Result<()> {
        let kmers = representation_methods::generate_kmers_with_step(
            &vec!['A', 'C', 'T', 'G'],
            2,
            2
        )?;
        let expected = vec![&['A', 'C'], &['T', 'G']];
        assert_eq!(kmers, expected);
        Ok(())
    }

    #[test]
    fn test_generate_g_spaced_kmers_with_step1() -> Result<()> {
        let gapmers = representation_methods::generate_g_spaced_kmers_with_step(
            &vec!['A', 'C', 'T', 'G', 'A', 'C', 'T'],
            3,
            1,
            1
        )?;
        let mut expected = HashMap::new();
        expected.insert(
            vec!['1', '0', '1', '1'],
            vec![vec!['A', 'T', 'G'], vec!['C', 'G', 'A'], vec!['T', 'A', 'C'], vec!['G', 'C', 'T']]
        );
        expected.insert(
            vec!['1', '1', '0', '1'],
            vec![vec!['A', 'C', 'G'], vec!['C', 'T', 'A'], vec!['T', 'G', 'C'], vec!['G', 'A', 'T']]
        );
        assert_eq!(gapmers, expected);

        Ok(())
    }

    #[test]
    fn test_generate_g_spaced_kmers_with_step2() -> Result<()> {
        let gapmers = representation_methods::generate_g_spaced_kmers_with_step(
            &vec!['A', 'C', 'T', 'G', 'A', 'C', 'T', 'G'],
            3,
            2,
            1
        )?;
        let mut expected = HashMap::new();
        expected.insert(
            vec!['1', '0', '0', '1', '1'],
            vec![vec!['A', 'G', 'A'], vec!['C', 'A', 'C'], vec!['T', 'C', 'T'], vec!['G', 'T', 'G']]
        );
        expected.insert(
            vec!['1', '0', '1', '0', '1'],
            vec![vec!['A', 'T', 'A'], vec!['C', 'G', 'C'], vec!['T', 'A', 'T'], vec!['G', 'C', 'G']]
        );

        expected.insert(
            vec!['1', '1', '0', '0', '1'],
            vec![vec!['A', 'C', 'A'], vec!['C', 'T', 'C'], vec!['T', 'G', 'T'], vec!['G', 'A', 'G']]
        );

        assert_eq!(expected, gapmers);
        Ok(())
    }

}
