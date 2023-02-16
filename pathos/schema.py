schema_string = """type: object

required:
    - star
    - pricefilter
    - bowtie2
    - ncbi
    - threads
    - assembly
    - cdhitdup
    - lzwfilter
    - blastn
    - param_file
    - ete2_db
properties:
  star:
    type: object
    properties:
      starDB:
        type: string
        pattern: ^(/[^/]+)+/?$
      skip:
        type: boolean
  pricefilter:
    type: object
    properties:
      highQualPercent:
        type: integer
        minimum: 1
        maximum: 100
      calledPercent:
        type: integer
        minimum: 1
        maximum: 100
      highQualMin:
        type: number
        minimum: 0
        maximum: 1
    required:
      - highQualPercent
      - highQualMin
      - calledPercent
    additionalProperties: false
  bowtie2:
    type: object
    properties:
      bowtieDB:
        type: string
        pattern: ^(/[^/]+)+/?$
    additionalProperties: false
    required:
      - bowtieDB
  threads:
    type: integer
    minimum: 1
    maximum: 64
  assembly:
    type: object
    properties:
      assembler:
        enum:
          - abyss
          - ray2
      options:
        type: string
      kmer:
        type: integer
      minimum_contig_length:
          type: integer
          minimum: 1
          maximum: 10000000
    required:
      - assembler
      - options
      - minimum_contig_length
      - kmer
    additionalProperties: false
  cdhitdup:
    type: object
    properties:
      minDifference:
        type: integer
        minimum: 1
        maximum: 100
    additionalProperties: false
    required:
      - minDifference
  lzwfilter:
    type: object
    properties:
      maxCompressionScore:
        type: number
        minimum: 0
        maximum: 1
  ncbi:
    type: object
    properties:
      ntDB:
        type: string
        pattern: ^(/[^/]+)+$
      nrDB:
        type: string
        pattern: ^(/[^/]+)+$
      ktTaxonomy:
        type: string
        pattern: ^(/[^/]+)+/?$
    required:
      - ntDB
      - ktTaxonomy
  blastn:
    type: object
    properties:
      max_target_seqs:
        type: number
        minimum: 1
        maximum: 100
  param_file:
    type: string
    pattern: ^(/[^/]+)+$
  ete2_db:
    type: string
    pattern: ^(/[^/]+)+$"""
import yaml
schema = yaml.load(schema_string)
