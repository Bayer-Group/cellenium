resource "aws_ecr_repository" "study_import_batch_repository" {
  name                 = "${data.aws_caller_identity.current.account_id}-${var.ecr_repository_name}"
  image_tag_mutability = "MUTABLE"

  image_scanning_configuration {
    scan_on_push = true
  }
}


# build docker image
resource "docker_image" "study_import_batch_image" {
  name = "${aws_ecr_repository.study_import_batch_repository.repository_url}:latest"
  build {
    context    = abspath("${abspath(path.module)}/../../../")
    dockerfile = "./terraform/modules/batch/Dockerfile" # path relative to context
  }
  platform = "linux/amd64"
}

# push image to ecr repo
resource "docker_registry_image" "media-handler" {
  name = docker_image.study_import_batch_image.name

  depends_on = [docker_image.study_import_batch_image]
}


resource "aws_security_group" "cellenium_study_import_batch_security_group" {
  name        = "${data.aws_caller_identity.current.account_id}-${var.security_group_name}"
  description = "Security group for cellenium study import with AWS Batch"
  vpc_id      = var.vpc_id

  egress {
    from_port       = var.ec2_security_group_port
    to_port         = var.ec2_security_group_port
    protocol        = "TCP"
    security_groups = [var.ec2_security_group_id]
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
  compute_environment_name = "${data.aws_caller_identity.current.account_id}-${var.security_group_name}"

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
  name = "${data.aws_caller_identity.current.account_id}-${var.job_queue_name}"

  priority = 1
  state    = "ENABLED"

  compute_environments = [aws_batch_compute_environment.cellenium_study_import_compute_environment.arn]
}

resource "aws_batch_job_definition" "cellenium_study_import_job_definition" {
  name = "${data.aws_caller_identity.current.account_id}-${var.job_definition_name}"
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
        value = var.db_secret_name
      },
      {
        name  = "AWS"
        value = "true"
      }
    ]
  })
}


resource "aws_cloudwatch_event_rule" "failed_job_cloudwatch_event_rule" {
  name        = "${data.aws_caller_identity.current.account_id}-${var.failed_batch_cloudwatch_rule_name}"
  description = "calls a lambda function is a batch import job failed."


  event_pattern = jsonencode({
    "detail-type" : [
      "Batch Job State Change"
    ],
    "source" : [
      "aws.batch"
    ],
    "detail" : {
      "status" : [
        "FAILED"
      ]
    }
  })
}


resource "aws_cloudwatch_event_target" "cellenium_batch_job_failed_lambda_target" {
  rule      = aws_cloudwatch_event_rule.failed_job_cloudwatch_event_rule.name
  target_id = aws_cloudwatch_event_rule.failed_job_cloudwatch_event_rule.name
  arn       = var.failed_lambda_arn
}