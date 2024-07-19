extern crate alignment_free_methods;

#[cfg(test)]
mod tests {
    use pretty_assertions::assert_eq;
    use anyhow::Result;

    use alignment_free_methods::utils::representation_methods;

    #[test]
    fn test_generate_kmers() -> Result<()> {
        let k2mers = representation_methods::generate_kmers(
            &vec!['A', 'C', 'T', 'G'],
            2
        )?;
        let expected = vec![&['A', 'C'], &['C', 'T'], &['T', 'G']];
        assert_eq!(k2mers, expected);
        Ok(())
    }

    #[test]
    fn test_generate_kmers_with_step() -> Result<()> {
        let k2mers_step2 = representation_methods::generate_kmers_with_step(
            &vec!['A', 'C', 'T', 'G'],
            2,
            2
        )?;
        let expected = vec![&['A', 'C'], &['T', 'G']];
        assert_eq!(k2mers_step2, expected);
        Ok(())
    }

}
