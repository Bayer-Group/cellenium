resource "aws_ecr_repository" "study_import_batch_repository" {
  name                 = "${data.aws_caller_identity.current.account_id}-${var.ecr_repository_name}"
  image_tag_mutability = "MUTABLE"

  image_scanning_configuration {
    scan_on_push = true
  }
}


resource "null_resource" "study_import_batch_docker_image" {
  triggers = {
    file_hash_postgres_utils            = filebase64sha256("${abspath(path.module)}/../../../data_import/postgres_utils.py")
    file_hash_study_import              = filebase64sha256("${abspath(path.module)}/../../../data_import/study_import.py")
    file_hash_study_import_requirements = filebase64sha256("${abspath(path.module)}/../../../data_import/study_import_requirements.txt")
    file_hash_dockerfile                = filebase64sha256("${abspath(path.module)}/Dockerfile")
  }

  provisioner "local-exec" {
    command = "docker login --password ${data.aws_ecr_authorization_token.token.password} --username ${data.aws_ecr_authorization_token.token.user_name} ${data.aws_ecr_authorization_token.token.proxy_endpoint} && docker buildx build --platform linux/amd64 -t ${aws_ecr_repository.study_import_batch_repository.repository_url}:latest --push -f ${abspath(path.module)}/Dockerfile ${abspath(path.module)}/../../../"
  }

}


## build docker image
#resource "docker_image" "study_import_batch_image" {
#  name = "${aws_ecr_repository.study_import_batch_repository.repository_url}:latest"
#
#  triggers = {
#    file_hash_postgres_utils            = filebase64sha256("${abspath(path.module)}/../../../data_import/postgres_utils.py")
#    file_hash_study_import              = filebase64sha256("${abspath(path.module)}/../../../data_import/study_import.py")
#    file_hash_study_import_requirements = filebase64sha256("${abspath(path.module)}/../../../data_import/study_import_requirements.txt")
#    file_hash_dockerfile                = filebase64sha256("${abspath(path.module)}/Dockerfile")
#  }
#
#  build {
#    context    = abspath("${abspath(path.module)}/../../../")
#    dockerfile = "./terraform/modules/batch/Dockerfile" # path relative to context
#    platform   = "linux/amd64"
#  }
#  platform = "linux/amd64"
#}

## push image to ecr repo
#resource "docker_registry_image" "media-handler" {
#  name = "${aws_ecr_repository.study_import_batch_repository.repository_url}:latest"
#
#  depends_on = [docker_image.study_import_batch_image]
#}


resource "aws_security_group" "cellenium_study_import_batch_security_group" {
  name        = "${data.aws_caller_identity.current.account_id}-${var.security_group_name}"
  description = "Security group for cellenium study import with AWS Batch"
  vpc_id      = var.vpc_id
}


resource "aws_security_group_rule" "cellenium_study_import_batch_egress_security_group_rule" {
  type                     = "egress"
  from_port                = var.ec2_security_group_port
  to_port                  = var.ec2_security_group_port
  protocol                 = "TCP"
  source_security_group_id = var.ec2_security_group_id
  security_group_id        = aws_security_group.cellenium_study_import_batch_security_group.id
}


resource "aws_security_group_rule" "cellenium_study_import_batch_security_group_rule" {
  type              = "egress"
  from_port         = 0
  to_port           = 65535
  protocol          = "TCP"
  cidr_blocks       = ["0.0.0.0/0"]
  security_group_id = aws_security_group.cellenium_study_import_batch_security_group.id
}

resource "aws_security_group_rule" "cellenium_study_import_security_group_rule_ec2" {
  type                     = "ingress"
  from_port                = var.ec2_security_group_port
  to_port                  = var.ec2_security_group_port
  protocol                 = "TCP"
  source_security_group_id = aws_security_group.cellenium_study_import_batch_security_group.id
  security_group_id        = var.ec2_security_group_id
  description              = "${data.aws_caller_identity.current.account_id}-cellenium-batch-access"
}

#resource "aws_security_group_rule" "cellenium_study_import_batch_ingress_security_group_rule" {
#  type              = "ingress"
#  from_port         = 0
#  to_port           = 65535
#  protocol          = "TCP"
#  cidr_blocks       = ["0.0.0.0/0"]
#  security_group_id = aws_security_group.cellenium_study_import_batch_security_group.id
#}


resource "aws_batch_compute_environment" "cellenium_study_import_compute_environment" {
  compute_environment_name = "${data.aws_caller_identity.current.account_id}-${var.compute_environment_name}"

  compute_resources {
    type               = "FARGATE"
    max_vcpus          = 120
    security_group_ids = [aws_security_group.cellenium_study_import_batch_security_group.id]
    subnets            = var.subnet_ids
  }
  type         = "MANAGED"
  service_role = aws_iam_role.batch_execution_role.arn
}

resource "aws_batch_job_queue" "cellenium_study_import_job_queue" {
  name = "${data.aws_caller_identity.current.account_id}-${var.job_queue_name}"

  priority = 1
  state    = "ENABLED"

  compute_environments = [aws_batch_compute_environment.cellenium_study_import_compute_environment.arn]
  depends_on           = [aws_batch_compute_environment.cellenium_study_import_compute_environment]
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
    image      = "${aws_ecr_repository.study_import_batch_repository.repository_url}:latest"
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

    #    networkConfiguration = {
    #      assignPublicIp = "ENABLED"
    #    },

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
        value = var.db_secret_id
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