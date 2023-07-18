resource "aws_iam_role" "batch_execution_role" {
  name = "cellenium_study_import_batch_role"

  assume_role_policy = jsonencode({
    Version   = "2012-10-17"
    Statement = [
      {
        Effect    = "Allow",
        Principal = {
          Service = "batch.amazonaws.com"
        },
        Action = "sts:AssumeRole"
      },
      {
        Effect    = "Allow",
        Principal = {
          Service = "ecs-tasks.amazonaws.com"
        },
        Action = "sts:AssumeRole"
      }
    ]
  })

  inline_policy {
    name = "cellenium_study_import_batch_policy"

    policy = jsonencode({
      Version   = "2012-10-17"
      Statement = [
        {
          "Effect" = "Allow",
          "Action" = [
            "s3:GetObject",
            "s3:ListBucket"
          ],
          "Resource" = [
            aws_s3_bucket.cellenium_s3_bucket.arn,
          ]
        },
        {
          "Effect" = "Allow",
          "Action" = [
            "secretsmanager:GetSecretValue",
            "secretsmanager:DescribeSecret",
            "secretsmanager:ListSecretVersionIds"
          ],
          "Resource" = aws_secretsmanager_secret.db_secret.arn
        },
        {
          "Effect" = "Allow",
          "Action" = [
            "ecs:DeleteCluster", "ecs:DescribeClusters", "ecs:ListClusters", "ecs:DeleteCluster"
          ],
          "Resource" = ["*"]
        },
                {
          "Effect" = "Allow",
          "Action" = [
            "logs:DescribeLogGroups"
          ],
          "Resource" = ["*"]
        }
      ]
    })
  }
}

resource "aws_iam_role_policy_attachment" "batch_role_policy_attachment_ecr" {
  role       = aws_iam_role.batch_execution_role.name
  policy_arn = "arn:aws:iam::aws:policy/AmazonEC2ContainerRegistryReadOnly"
}

resource "aws_iam_role_policy_attachment" "batch_role_policy_attachment_batch_sr" {
  role       = aws_iam_role.batch_execution_role.name
  policy_arn = "arn:aws:iam::aws:policy/service-role/AWSBatchServiceRole"
}


resource "aws_iam_role_policy_attachment" "batch_role_policy_attachment_ecs_task" {
  role       = aws_iam_role.batch_execution_role.name
  policy_arn = "arn:aws:iam::aws:policy/service-role/AmazonECSTaskExecutionRolePolicy"
}


resource "aws_ecr_repository" "study_import_batch_repository" {
  name                 = "cellenium_study_import_batch"
  image_tag_mutability = "MUTABLE"

  image_scanning_configuration {
    scan_on_push = true
  }
}


# build docker image
resource "docker_image" "study_import_batch_image" {
  name = "${aws_ecr_repository.study_import_batch_repository.repository_url}:latest"
  build {
    context    = abspath("${abspath(path.cwd)}/../")
    dockerfile = "pipeline/batch/Dockerfile"
  }
  platform = "linux/amd64"
}

# push image to ecr repo
resource "docker_registry_image" "media-handler" {
  name = docker_image.study_import_batch_image.name

  depends_on = [docker_image.study_import_batch_image]
}


resource "aws_security_group" "cellenium_study_import_batch_security_group" {
  name        = "cellenium_study_import_batch_security_group"
  description = "Security group for cellenium study import with AWS Batch"
  vpc_id      = data.aws_vpc.vpc.id

  egress {
    from_port       = var.cellenium_ec2_access_security_group_port
    to_port         = var.cellenium_ec2_access_security_group_port
    protocol        = "TCP"
    security_groups = [var.cellenium_ec2_access_security_group_id]
  }
}

resource "aws_security_group_rule" "cellenium_study_import_batch_security_group_rule" {
  type              = "egress"
  from_port         = 0
  to_port           = 65535
  protocol          = "TCP"
  cidr_blocks       = [data.aws_vpc.vpc.cidr_block]
  security_group_id = aws_security_group.cellenium_study_import_batch_security_group.id
}


resource "aws_batch_compute_environment" "cellenium_study_import_compute_environment" {
  compute_environment_name = "cellenium_study_import_compute_environment"

  compute_resources {
    type               = "FARGATE"
    max_vcpus          = 120
    security_group_ids = [aws_security_group.cellenium_study_import_batch_security_group.id]
    subnets            = data.aws_subnets.default_vpc_subnets.ids
  }
  type         = "MANAGED"
  service_role = aws_iam_role.batch_execution_role.arn
}

resource "aws_batch_job_queue" "cellenium_study_import_job_queue" {
  name = "cellenium_study_import_job_queue"

  priority = 1
  state    = "ENABLED"

  compute_environments = [aws_batch_compute_environment.cellenium_study_import_compute_environment.arn]
}

resource "aws_batch_job_definition" "cellenium_study_import_job_definition" {
  name = "cellenium_study_import"
  type = "container"

  platform_capabilities = [
    "FARGATE",
  ]

  timeout {
    attempt_duration_seconds = 86400
  }

  container_properties = jsonencode({

    command    = ["Ref::analyze-database", "Ref::filename"]
    image      = docker_image.study_import_batch_image.name
    jobRoleArn = aws_iam_role.batch_execution_role.arn

    fargatePlatformConfiguration = {
      platformVersion = "LATEST"
    }

    parameters = {
      analyze-database = "--analyze-database"
    }

    executionRoleArn = aws_iam_role.batch_execution_role.arn

    "ephemeralStorage" : {
      "sizeInGiB" : 100
    },

    networkConfiguration = {
      assignPublicIp = "ENABLED"
    },

    resourceRequirements = [
      {
        type  = "VCPU"
        value = "16"
      },
      {
        type  = "MEMORY"
        value = "65536"
      }
    ]

    environment = [
      {
        name  = "AWS_DB_SECRET"
        value = var.cellenium_db_secret_name
      },
      {
        name  = "AWS"
        value = "true"
      }
    ]
  })
}

