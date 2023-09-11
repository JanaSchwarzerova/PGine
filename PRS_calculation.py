"""
PRS_calculation.py - this script shows a short example of Polygenic Risk Score (PRS) calculation
This PRS calculation is applied in the main script to real plant data (Maize)
"""

# Sample genetic data (variant alleles)
genetic_data = {
    'rs12345': 2,  # Number of risk alleles (0, 1, or 2)
    'rs56789': 1,
    'rs98765': 0,
    # Add more variants as needed
}

# Effect sizes (beta coefficients) for each variant
effect_sizes = {
    'rs12345': 0.1,
    'rs56789': 0.09,
    'rs98765': -0.75,
    # Match each variant in genetic_data with its effect size
}

# Calculate the Polygenic Risk Score
prs = 0
for variant, alleles in genetic_data.items():
    if variant in effect_sizes:
        prs += effect_sizes[variant] * alleles

print("Polygenic Risk Score (PRS):", prs)
