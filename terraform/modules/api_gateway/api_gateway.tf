resource "aws_lb" "load_balancer" {
  name               = "${data.aws_caller_identity.current.account_id}-${var.stage}-${var.nlb_name}"
  load_balancer_type = "network"
  subnets            = var.nlb_subnet_ids
  internal           = true

  enable_deletion_protection = true
  security_groups            = []
}


resource "aws_lb_listener" "nlb_listener" {
  load_balancer_arn = aws_lb.load_balancer.arn
  port              = var.nlb_port
  protocol          = "TCP"

  default_action {
    type             = "forward"
    target_group_arn = aws_lb_target_group.nlb_target_group.arn
  }
}

resource "aws_security_group" "cellenium_nlb_security_group" {
  name   = "${data.aws_caller_identity.current.account_id}-${var.stage}-${var.nlb_security_group_name}"
  vpc_id = var.vpc_id

  egress {
    from_port       = 0
    to_port         = 0
    protocol        = "-1"
    security_groups = [var.ec2_security_group_id]
  }


  ingress {
    from_port       = 0
    to_port         = 0
    protocol        = "-1"
    security_groups = [aws_security_group.cellenium_vpc_link_security_group.id]
  }
}


resource "aws_lb_target_group" "nlb_target_group" {
  name             = "${var.stage}-${var.nlb_target_group_name}"
  port             = var.ec2_target_group_port
  protocol         = "HTTP"
  protocol_version = "HTTP1"
  vpc_id           = var.vpc_id
  target_type      = "instance"
}


resource "aws_lb_target_group_attachment" "test" {
  target_group_arn = aws_lb_target_group.nlb_target_group.arn
  target_id        = var.ec2_instance_id
  port             = var.ec2_target_group_port
}


resource "aws_security_group_rule" "ec2_security_group_rule" {
  type                     = "ingress"
  from_port                = var.ec2_security_group_port
  to_port                  = var.ec2_security_group_port
  protocol                 = "TCP"
  source_security_group_id = aws_security_group.cellenium_nlb_security_group.id
  security_group_id        = var.ec2_security_group_id
}


resource "aws_apigatewayv2_api" "api_gateway" {
  name          = "${data.aws_caller_identity.current.account_id}-${var.stage}-${var.api_name}"
  protocol_type = "HTTP"
}

resource "aws_cloudwatch_log_group" "cloudwatch_log_group" {
  name = "${data.aws_caller_identity.current.account_id}-${var.stage}-${var.api_log_group_name}"
}

resource "aws_apigatewayv2_stage" "example" {
  api_id = aws_apigatewayv2_api.api_gateway.id
  name   = "default"

  access_log_settings {
    destination_arn = aws_cloudwatch_log_group.cloudwatch_log_group.arn
    format          = "$context.identity.sourceIp,$context.requestTime,$context.httpMethod,$context.path,$context.routeKey,$context.protocol,$context.status,$context.responseLength,$context.requestId"
  }

  default_route_settings {
    logging_level            = "INFO"
    detailed_metrics_enabled = true
  }
}

resource "aws_apigatewayv2_authorizer" "oidc_authorizer" {
  api_id           = aws_apigatewayv2_api.api_gateway.id
  authorizer_type  = "JWT"
  identity_sources = ["$request.header.Authorization"]
  name             = "${data.aws_caller_identity.current.account_id}-${var.stage}-${var.api_authorizer_name}"

  jwt_configuration {
    audience = var.api_authorizer_audiences
    issuer   = var.api_authorizer_issuer
  }
}

resource "aws_apigatewayv2_domain_name" "api" {
  domain_name = var.api_domain_name

  domain_name_configuration {
    certificate_arn = var.api_domain_name_certificate_arn
    endpoint_type   = "REGIONAL"
    security_policy = "TLS_1_2"
  }
}

resource "aws_route53_record" "example" {
  name    = var.api_domain_name
  type    = "A"
  zone_id = var.route53_hosted_zone_id

  alias {
    name                   = aws_apigatewayv2_domain_name.api.domain_name_configuration[0].target_domain_name
    zone_id                = aws_apigatewayv2_domain_name.api.domain_name_configuration[0].hosted_zone_id
    evaluate_target_health = false
  }
}

resource "aws_apigatewayv2_vpc_link" "vpc_link" {
  name               = "${data.aws_caller_identity.current.account_id}-${var.stage}-${var.vpc_link_name}"
  security_group_ids = [aws_security_group.cellenium_vpc_link_security_group.id]
  subnet_ids         = var.vpc_link_subnet_ids
}


resource "aws_security_group" "cellenium_vpc_link_security_group" {
  name   = "${data.aws_caller_identity.current.account_id}-${var.stage}-${var.vpc_link_security_group_name}"
  vpc_id = var.vpc_id

  egress {
    from_port        = 0
    to_port          = 0
    protocol         = "-1"
    cidr_blocks      = ["0.0.0.0/0"]
    ipv6_cidr_blocks = ["::/0"]
  }
}


resource "aws_apigatewayv2_route" "example" {
  api_id    = aws_apigatewayv2_api.api_gateway.id
  route_key = "$default"

}

resource "aws_apigatewayv2_integration" "example" {
  api_id           = aws_apigatewayv2_api.api_gateway.id
  description      = "cellenium api load balancer integration"
  integration_type = "HTTP_PROXY"
  integration_uri  = aws_lb_listener.nlb_listener.arn

  integration_method = "ANY"
  connection_type    = "VPC_LINK"
  connection_id      = aws_apigatewayv2_vpc_link.vpc_link.id

  tls_config {
    server_name_to_verify = var.api_domain_name
  }

  request_parameters = {
    "append:header.authforintegration" = "$context.authorizer.authorizerResponse"
    "overwrite:path"                   = "staticValueForIntegration"
  }

  response_parameters {
    status_code = 403
    mappings    = {
      "append:header.auth" = "$context.authorizer.authorizerResponse"
    }
  }
}